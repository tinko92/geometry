#!/usr/bin/env python3
"""Exclusive Callgrind Ir analysis; no profiler, compiler, or target is launched.

python3 analyze_simple_cartesian_profile.py DATA --medians native-analysis.json \
    --output /tmp/profile-analysis.json --text /tmp/profile-analysis.txt
The parser follows https://valgrind.org/docs/manual/cl-format.html.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import statistics


GROUPS = ("predicate_filter", "construction_filter", "checked_backend_validation", "exact_arithmetic",
          "numerical_shared", "representation_support", "topology_parenting",
          "memory", "allocation", "bgl", "algorithm", "standard_library",
          "harness", "unknown_system", "unknown")
ALLOC = re.compile(r"(?:^|::|\s)(?:malloc|calloc|realloc|free|cfree|aligned_alloc|memalign|posix_memalign|_int_malloc|_int_free)(?:$|\s)|operator (?:new|delete)")
MEMORY = re.compile(r"mem(?:cpy|move|set|cmp)|bzero|(?:^|::)(?:__copy|__fill|copy|copy_n|fill|fill_n|uninitialized_copy|uninitialized_move)(?:_|$|::)")


def stem(symbol):
    """Remove template arguments and function arguments, not just their text."""
    if symbol.startswith(("decltype(", "__decltype(")):
        depth = 0
        for i in range(symbol.index("("), len(symbol)):
            depth += (symbol[i] == "(") - (symbol[i] == ")")
            if depth == 0:
                symbol = symbol[i + 1:].lstrip()
                break
    out, depth = [], 0
    for c in symbol:
        if c == "<":
            depth += 1
        elif c == ">" and depth:
            depth -= 1
        elif not depth:
            if c == "(":
                break
            out.append(c)
    return " ".join("".join(out).split())


def short_file(filename):
    marker = "/include/boost/"
    if marker in filename:
        return "boost/" + filename.split(marker, 1)[1]
    return filename


def unresolved(name):
    return name == "???" or bool(re.fullmatch(r"0x[0-9a-fA-F]+", name))


def classify(filename, name, obj, allocator_inferred=False):
    """Use the actual line's file and the argument-free function owner/name."""
    base = Path(filename).name
    if ALLOC.search(name) or allocator_inferred or base in ("new_allocator.h", "malloc.c"):
        return "allocation"
    if MEMORY.search(name) or base.startswith(("memcpy", "memmove", "memset")):
        return "memory"
    if "multiprecision" in filename and re.search(r"::verify_(?:new_size|limb_mask)$", name):
        return "checked_backend_validation"
    if base in ("fixed_arithmetic.hpp", "fixed_rounding.hpp", "fixed_measures.hpp",
                "construction_filter.hpp", "static_number.hpp") or "multiprecision" in filename:
        return "exact_arithmetic"
    if base == "construction_rounding.hpp" or "::native_intersection" in name or "::segment_parameter" in name:
        return "construction_filter"
    if base in ("borrowed_predicates.hpp", "ulp_filter.hpp", "static_predicates.hpp",
                "orientation.hpp", "side_robust.hpp", "precise_math.hpp", "native_predicates.hpp",
                "segment_position.hpp", "static_box.hpp"):
        return "predicate_filter"
    if base == "forward_error.hpp":
        return "predicate_filter" if "::try_" in name else "numerical_shared"
    if base in ("arithmetic.hpp", "two_product.hpp", "simple_cartesian.hpp", "cmath") or name.startswith(("fma", "nextafter")):
        return "numerical_shared"
    if base in ("ring_topology.hpp", "set_ring_parents.hpp", "ring_properties.hpp",
                "topological_ring.hpp", "static_ring_topology.hpp", "borrowed_ring.hpp") or "::ring_topology_detail::" in name:
        return "topology_parenting"
    if base in ("borrowed_records.hpp", "borrowed_bounds.hpp", "borrowed_geometry.hpp",
                "borrowed_strategy.hpp", "point_cache.hpp", "operation_context.hpp",
                "query_point.hpp", "point.hpp", "borrowed_query.hpp", "borrowed_native.hpp",
                "borrowed_output.hpp", "borrowed_turns.hpp", "iterator_point_view.hpp",
                "point_precision_cache.hpp", "native_traits.hpp", "point_xy.hpp", "variant"):
        return "representation_support"
    if "boost/graph/" in filename or "::adjacency_list" in name or "::disjoint_sets" in name:
        return "bgl"
    if "boost/geometry/" in filename or "boost::geometry::" in name:
        return "algorithm"
    if "/c++/" in filename or name.startswith(("std::", "void std::", "bool std::")):
        return "standard_library"
    if base in ("profile_simple_cartesian_calls.cpp", "simple_cartesian_benchmark.cpp") or name in ("main", "execute"):
        return "harness"
    if unresolved(name) and (".so" in obj or filename == "???"):
        return "unknown_system"
    return "unknown"


class Callgrind:
    def __init__(self):
        self.header = {}
        self.names = {"file": {}, "fn": {}, "object": {}}
        self.functions = []
        self.function_ids = {}
        self.lines = Counter()
        self.arcs = Counter()
        self.arc_ir = Counter()
        self.summary = None
        self.footer = None

    def name(self, family, value):
        match = re.fullmatch(r"\((\d+)\)(?:\s+(.*))?", value)
        if not match:
            return value
        number, text = match.groups()
        if text is not None:
            self.names[family][number] = text
        if number not in self.names[family]:
            raise ValueError("Undefined compressed name: " + value)
        return self.names[family][number]

    def function(self, obj, symbol, filename):
        key = (obj, symbol)
        if key not in self.function_ids:
            self.function_ids[key] = len(self.functions)
            self.functions.append({"object": obj, "name": stem(symbol),
                                   "symbol_sha256": hashlib.sha256(symbol.encode()).hexdigest(),
                                   "definition_file": short_file(filename)})
        return self.function_ids[key]

    def read(self, text):
        obj = filename = symbol = "???"
        callee_obj = callee_file = callee_symbol = "???"
        positions, events, previous, pending, function = ["line"], [], [0], None, None
        parts = 0
        for lineno, raw in enumerate(text.splitlines(), 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if re.match(r"^[A-Za-z_]+:", line):
                key, value = line.split(":", 1)
                value = value.strip()
                if key == "part":
                    parts += 1
                    if parts > 1:
                        raise ValueError("Multiple parts require separate input files")
                if key == "positions":
                    positions, previous = value.split(), [0] * len(value.split())
                elif key == "events":
                    events = value.split()
                    if "Ir" not in events:
                        raise ValueError("Profile does not record Ir")
                elif key in ("summary", "totals"):
                    values = value.split()
                    count = int(values[events.index("Ir")])
                    if key == "summary": self.summary = count
                    else: self.footer = count
                self.header[key] = value
                continue
            if "=" in line:
                key, value = line.split("=", 1)
                value = value.strip()
                if key in ("ob", "cob"):
                    resolved = self.name("object", value)
                    if key == "ob": obj = resolved
                    else: callee_obj = resolved
                elif key in ("fl", "fi", "fe", "cfi", "cfl", "jfi"):
                    resolved = self.name("file", value)
                    if key in ("fl", "fi", "fe"): filename = resolved
                    elif key in ("cfi", "cfl"): callee_file = resolved
                elif key in ("fn", "cfn"):
                    resolved = self.name("fn", value)
                    if key == "fn":
                        symbol = resolved
                        function = self.function(obj, symbol, filename)
                        callee_obj, callee_file, callee_symbol = obj, filename, "???"
                    else: callee_symbol = resolved
                elif key == "calls":
                    if pending is not None:
                        raise ValueError("Call association lacks cost line")
                    pending = (int(value.split()[0]), self.function(callee_obj, callee_symbol, callee_file))
                elif key not in ("jump", "jcnd", "jfl", "jfn"):
                    raise ValueError("Unknown position record: " + key)
                continue
            fields = line.split()
            if function is None or not events or len(fields) < len(positions):
                raise ValueError("Invalid cost line at " + str(lineno))
            for i, token in enumerate(fields[:len(positions)]):
                if token != "*":
                    value = int(token.lstrip("+-"), 16 if token.lstrip("+-").startswith("0x") else 10)
                    previous[i] = previous[i] + (-value if token.startswith("-") else value) if token[0] in "+-" else value
            costs = fields[len(positions):]
            index = events.index("Ir")
            ir = int(costs[index]) if index < len(costs) else 0
            source_line = previous[positions.index("line")] if "line" in positions else 0
            if pending is not None:
                count, target = pending
                arc = (function, target, short_file(filename), source_line)
                self.arcs[arc] += count
                self.arc_ir[arc] += ir
                pending = None
            else:
                self.lines[(function, short_file(filename), source_line)] += ir
        if pending is not None:
            raise ValueError("Unterminated call association")
        self.total = sum(self.lines.values())
        if self.footer is None or self.total != self.footer:
            raise ValueError("Exclusive Ir sum does not equal totals footer: %s != %s" % (self.total, self.footer))
        return self

    def allocator_closure(self):
        incoming = defaultdict(set)
        for (caller, callee, _, _), calls in self.arcs.items():
            if calls: incoming[callee].add(caller)
        inferred = {}
        for target, fn in enumerate(self.functions):
            if not unresolved(fn["name"]): continue
            pending, visited, frontier = list(incoming[target]), {target}, set()
            while pending:
                caller = pending.pop()
                if caller in visited: continue
                visited.add(caller)
                if unresolved(self.functions[caller]["name"]): pending.extend(incoming[caller])
                else: frontier.add(caller)
            if frontier and all(ALLOC.search(self.functions[i]["name"]) for i in frontier):
                inferred[target] = sorted(frontier)
        return inferred


def call_family(name):
    if ALLOC.search(name): return "allocator_api_layers_not_additive"
    if "refine" in name or "precision_homogeneous" in name: return "refinement"
    if re.search(r"same_intersection|same_.*support|on_support|shared_support|producing_support", name): return "shared_support"
    if "verify_new_size" in name or "verify_limb_mask" in name: return "checked_backend_validation"
    if "::try_compare" in name or "::try_equal" in name: return "point_filters"
    if "orientation" in name or "orient2d" in name: return "orientation"
    if "construction_rounding" in name or "native_intersection" in name: return "construction"
    if "::producers" in name or "::native_coordinate" in name: return "producer_access"
    return None


def analyze_profile(file, medians):
    profile = Callgrind().read(file.read_text())
    command = shlex.split(profile.header.get("cmd", ""))
    if len(command) < 4: raise ValueError("Missing wrapper command arguments")
    variant, category, operation, rounds = Path(command[-4]).name, command[-3], command[-2], int(command[-1])
    calls = rounds * 3
    if rounds <= 0: raise ValueError("Nonpositive rounds")
    inferred = profile.allocator_closure()
    group_costs, file_costs, function_costs, incoming_calls = Counter(), Counter(), Counter(), Counter()
    line_rows, arc_rows = [], []
    for (caller, callee, filename, line), count in profile.arcs.items():
        incoming_calls[callee] += count
        arc_rows.append({"caller": caller, "callee": callee, "file": filename, "line": line,
                         "calls": count, "calls_per_operation": count / calls,
                         "inclusive_ir_not_additive": profile.arc_ir[(caller, callee, filename, line)]})
    for (function, filename, line), ir in profile.lines.items():
        if not ir: continue
        fn = profile.functions[function]
        group = classify(filename, fn["name"], fn["object"], function in inferred)
        group_costs[group] += ir
        file_costs[filename] += ir
        function_costs[function] += ir
        line_rows.append({"function": function, "file": filename, "line": line,
                          "group": group, "self_ir": ir, "ir_per_call": ir / calls})
    assert sum(group_costs.values()) == profile.total
    for i, fn in enumerate(profile.functions):
        fn.update(id=i, self_ir=function_costs[i], self_ir_per_call=function_costs[i]/calls,
                  incoming_calls=incoming_calls[i], incoming_calls_per_operation=incoming_calls[i]/calls,
                  call_family=call_family(fn["name"]))
        if i in inferred:
            fn["allocation_inference"] = {"rule": "all resolved incoming caller frontiers are allocator APIs",
                                          "frontier_function_ids": inferred[i]}
    observed_execute_calls = sum(count for (_, callee, _, _), count in profile.arcs.items()
                                 if profile.functions[callee]["name"].split()[-1] == "execute")
    if observed_execute_calls != calls:
        raise ValueError("execute call count %d != wrapper denominator %d" % (observed_execute_calls, calls))
    annotation = file.with_name(file.name + "-functions.txt")
    annotated_total = None
    if annotation.exists():
        match = re.search(r"^\s*([\d,]+)\s+PROGRAM TOTALS\s*$", annotation.read_text(), re.M)
        if not match: raise ValueError("Missing annotated PROGRAM TOTALS")
        annotated_total = int(match[1].replace(",", ""))
        if annotated_total != profile.total: raise ValueError("Annotation total mismatch")
    median = medians.get((variant, category, operation))
    invalid = category in ("general", "decimal") and operation == "sym_difference"
    return {"file": str(file), "sha256": hashlib.sha256(file.read_bytes()).hexdigest(),
            "variant": variant, "category": category, "operation": operation, "rounds": rounds,
            "calls": calls, "execute_arc_calls": observed_execute_calls,
            "summary_header_ir": profile.summary, "totals_footer_ir": profile.footer,
            "annotated_total_ir": annotated_total, "self_ir": profile.total, "ir_per_call": profile.total/calls,
            "native_median_ns_per_call": median, "comparison_qualification":
            "not like-for-like: known invalid simple output; retained separately" if invalid else "controls diagnostic, not an oracle",
            "groups": {g: {"self_ir": group_costs[g], "ir_per_call": group_costs[g]/calls,
                            "instruction_share": group_costs[g]/profile.total} for g in GROUPS},
            "allocation_inferred_ir": sum(function_costs[i] for i in inferred),
            "top_files": [{"file": f, "self_ir": n, "ir_per_call": n/calls} for f, n in file_costs.most_common()],
            "functions": profile.functions, "self_lines": sorted(line_rows, key=lambda r: -r["self_ir"]),
            "call_arcs": sorted(arc_rows, key=lambda r: -r["calls"])}


def load_medians(file):
    if file is None: return {}
    data = json.loads(file.read_text())
    result = {}
    for row in data["cases"]:
        for side in ("candidate", "reference"):
            result[(row[side], row["category"], row["operation"])] = row[side + "_ns"]
    return result


def comparisons(profiles):
    by_key = {(p["variant"], p["category"], p["operation"]): p for p in profiles}
    result = []
    for (variant, category, operation), exact in sorted(by_key.items()):
        if variant != "exact" or ("simple", category, operation) not in by_key: continue
        simple = by_key[("simple", category, operation)]
        a, b = exact["native_median_ns_per_call"], simple["native_median_ns_per_call"]
        result.append({"category": category, "operation": operation,
                       "primary_25_case_comparison": not (category in ("general", "decimal") and operation == "sym_difference"),
                       "exact_ir_per_call": exact["ir_per_call"], "simple_ir_per_call": simple["ir_per_call"],
                       "ir_ratio": exact["ir_per_call"]/simple["ir_per_call"],
                       "exact_native_ns": a, "simple_native_ns": b,
                       "native_time_ratio": a/b if a is not None and b is not None else None,
                       "groups": {g: {"exact_ir_per_call": exact["groups"][g]["ir_per_call"],
                                      "simple_ir_per_call": simple["groups"][g]["ir_per_call"],
                                      "delta_ir_per_call": exact["groups"][g]["ir_per_call"]-simple["groups"][g]["ir_per_call"]}
                                  for g in GROUPS}})
    return result


def aggregate_cases(cases):
    if not cases: return {"cases": 0}
    count = len(cases)
    a = sum(r["exact_ir_per_call"] for r in cases) / count
    b = sum(r["simple_ir_per_call"] for r in cases) / count
    result = {"cases": count, "mean_exact_ir_per_call": a, "mean_simple_ir_per_call": b,
              "geomean_ir_ratio": math.exp(statistics.mean(math.log(r["ir_ratio"]) for r in cases)),
              "groups": {}}
    timed = [r for r in cases if r["native_time_ratio"] is not None]
    if timed:
        result.update(geomean_native_time_ratio=math.exp(statistics.mean(math.log(r["native_time_ratio"]) for r in timed)),
                      mean_exact_native_ns_per_call=statistics.mean(r["exact_native_ns"] for r in timed),
                      mean_simple_native_ns_per_call=statistics.mean(r["simple_native_ns"] for r in timed))
    for group in GROUPS:
        x = sum(r["groups"][group]["exact_ir_per_call"] for r in cases) / count
        y = sum(r["groups"][group]["simple_ir_per_call"] for r in cases) / count
        result["groups"][group] = {"exact_ir_per_call": x, "simple_ir_per_call": y,
                                    "delta_ir_per_call": x-y, "exact_instruction_share": x/a,
                                    "simple_instruction_share": y/b, "share_of_ir_delta": (x-y)/(a-b)}
    return result


def primary_hotspots(profiles):
    result = {}
    for variant in ("exact", "simple"):
        selected = [p for p in profiles if p["variant"] == variant
                    and not (p["operation"] == "sym_difference" and p["category"] in ("general", "decimal"))]
        if not selected: continue
        files, lines, functions, calls, families = Counter(), Counter(), Counter(), Counter(), defaultdict(Counter)
        inferred = 0
        for p in selected:
            inferred += p["allocation_inferred_ir"] / p["calls"]
            for r in p["top_files"]: files[r["file"]] += r["ir_per_call"]
            for r in p["self_lines"]: lines[(r["file"], r["line"], r["group"])] += r["ir_per_call"]
            for f in p["functions"]:
                functions[f["name"]] += f["self_ir_per_call"]
                calls[f["name"]] += f["incoming_calls_per_operation"]
                if f["call_family"]: families[f["call_family"]][f["name"]] += f["incoming_calls_per_operation"]
        n = len(selected)
        result[variant] = {
            "cases": n, "allocation_inferred_mean_ir_per_call": inferred/n,
            "top_files": [{"file": f, "mean_ir_per_call": ir/n} for f, ir in files.most_common(25)],
            "hot_lines": [{"file": f, "line": line, "group": g, "mean_ir_per_call": ir/n}
                          for (f, line, g), ir in lines.most_common(30)],
            "top_functions": [{"name": name, "mean_self_ir_per_call": ir/n,
                               "observed_calls_per_operation": calls[name]/n}
                              for name, ir in functions.most_common(25)],
            "observed_call_families_not_additive": {
                family: [{"name": name, "calls_per_operation": c/n} for name, c in counter.most_common(20)]
                for family, counter in families.items()}}
    return result


def rounded_operation_metrics(profiles):
    """Overlapping focus metric, not another disjoint attribution group."""
    filename = "boost/geometry/strategies/cartesian/exact/arithmetic.hpp"
    result = {"file": filename, "lines": {84: "rounded_add", 85: "rounded_sub", 86: "rounded_mul", 87: "rounded_div"},
              "meaning": "Self Ir attributed to entire source lines, including arithmetic/store/load/return; not isolated spill overhead or estimated time savings. Subset of numerical_shared.",
              "variants": {}}
    for variant in ("exact", "simple"):
        selected = [p for p in profiles if p["variant"] == variant and not
                    (p["operation"] == "sym_difference" and p["category"] in ("general", "decimal"))]
        if not selected: continue
        costs = Counter()
        for p in selected:
            for r in p["self_lines"]:
                if r["file"] == filename and r["line"] in result["lines"]:
                    costs[r["line"]] += r["ir_per_call"]
        n = len(selected)
        result["variants"][variant] = {"primary_cases": n,
            "mean_ir_per_call_by_line": {line: costs[line]/n for line in result["lines"]},
            "mean_ir_per_call": sum(costs.values())/n,
            "share_of_variant_primary_ir": sum(costs.values())/sum(p["ir_per_call"] for p in selected)}
    return result


def render(data):
    lines = [data["meaning"], "", "VALIDATION " + json.dumps(data["assurances"])]
    for key in ("primary_25", "invalid_simple_2", "all_27_retained"):
        summary = data["aggregates"][key]
        lines += ["", key + ": " + str(summary["cases"])]
        if not summary["cases"]: continue
        lines.append("Ir ratio GM %.4f; mean Ir/call exact %.2f simple %.2f; native time ratio GM %s" %
                     (summary["geomean_ir_ratio"], summary["mean_exact_ir_per_call"], summary["mean_simple_ir_per_call"],
                      summary.get("geomean_native_time_ratio")))
        for group in GROUPS:
            r = summary["groups"][group]
            lines.append("%s mean Ir/call exact=%.2f simple=%.2f delta=%.2f; exact share %.2f%%" %
                         (group, r["exact_ir_per_call"], r["simple_ir_per_call"], r["delta_ir_per_call"], 100*r["exact_instruction_share"]))
    for variant, summary in data["primary_hotspots"].items():
        lines += ["", variant + " PRIMARY TOP FILES"]
        lines += ["%.2f %s" % (r["mean_ir_per_call"], r["file"]) for r in summary["top_files"][:10]]
        lines += [variant + " PRIMARY HOT LINES"]
        lines += ["%.2f %s:%d %s" % (r["mean_ir_per_call"], r["file"], r["line"], r["group"]) for r in summary["hot_lines"][:10]]
        lines += [variant + " OBSERVED CALLS " + json.dumps(summary["observed_call_families_not_additive"])]
    lines += ["", "ROUNDED OPERATION LINES 84-87 " + json.dumps(data["rounded_operation_lines"])]
    lines += ["", "category operation primary Ir_ratio exact_Ir/call simple_Ir/call exact_ns simple_ns native_time_ratio"]
    for r in data["comparisons"]:
        lines.append("%s %s %s %.4f %.2f %.2f %s %s %s" % (r["category"], r["operation"], r["primary_25_case_comparison"],
                     r["ir_ratio"], r["exact_ir_per_call"], r["simple_ir_per_call"], r["exact_native_ns"], r["simple_native_ns"], r["native_time_ratio"]))
    return "\n".join(lines) + "\n"


def self_test():
    text = """events: Ir
positions: line
summary: 0
ob=(1) libc.so
fl=(1) ???
fn=(1) free
0 3
cfn=(2) 0x123
calls=2 0
* 7
fn=(2)
0 7
fl=(2) /include/boost/geometry/strategies/cartesian/exact/forward_error.hpp
fn=(3) bool boost::geometry::detail::compare<borrowed_root<double> >(int)
10 5
fi=(3) inline.hpp
+2 2
fe=(2)
* 1
jcnd=2 1 0
+1 4
totals: 22
"""
    p = Callgrind().read(text)
    assert p.total == 22 and p.summary == 0 and sum(p.arc_ir.values()) == 7
    assert p.allocator_closure() == {1: [0]}
    assert p.lines[(2, "inline.hpp", 12)] == 2
    assert "borrowed_root" not in stem("void get_turns<borrowed_root<double> >::apply(native_point)")
    assert stem("decltype(auto) ns::visit<Root>(Point)") == "ns::visit"
    assert classify("boost/geometry/algorithms/detail/overlay/get_turns.hpp", "void get_turns::apply", "a") == "algorithm"
    try: Callgrind().read(text.replace("totals: 22", "totals: 23"))
    except ValueError: pass
    else: raise AssertionError("Must reject inconsistent totals")
    print("Callgrind parser self-test passed.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path)
    parser.add_argument("--medians", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--text", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test: self_test(); return
    files = []
    for file in args.inputs:
        files.extend(sorted(file.glob("*.callgrind")) if file.is_dir() else [file])
    if not files: parser.error("Provide Callgrind file(s) or directory")
    medians = load_medians(args.medians)
    profiles = [analyze_profile(file, medians) for file in files]
    result = {"meaning": "Exclusive Ir is executed instruction work, NOT cycles or proportional CPU time. Native uninstrumented median ns/call is shown separately. Call-arc inclusive costs are not added to self costs. Calls count only observable out-of-line arcs; inlined predicates/refinements may have no separate call record.",
              "format_specification": "https://valgrind.org/docs/manual/cl-format.html",
              "assurances": {"profiles": len(profiles), "operation_calls": sum(p["calls"] for p in profiles),
                             "all_exclusive_sums_match_totals": True, "all_execute_call_denominators_match": True,
                             "annotated_totals_checked": sum(p["annotated_total_ir"] is not None for p in profiles)},
              "profiles": profiles, "comparisons": comparisons(profiles)}
    cases = result["comparisons"]
    result["aggregates"] = {"primary_25": aggregate_cases([r for r in cases if r["primary_25_case_comparison"]]),
                            "invalid_simple_2": aggregate_cases([r for r in cases if not r["primary_25_case_comparison"]]),
                            "all_27_retained": aggregate_cases(cases)}
    result["primary_hotspots"] = primary_hotspots(profiles)
    result["rounded_operation_lines"] = rounded_operation_metrics(profiles)
    source = Path(__file__).resolve().parents[2] / "include/boost/geometry/strategies/cartesian/exact/arithmetic.hpp"
    if source.exists():
        result["rounded_operation_lines"]["local_source_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    result["timing_reference"] = {"file": str(args.medians), "sha256": hashlib.sha256(args.medians.read_bytes()).hexdigest()} if args.medians else None
    result["classification_notes"] = [
        "Groups use exclusive inline source-file costs and argument-free owner/function stems; Root template arguments never classify a sample.",
        "Plain simple_cartesian strategy code is numerical_shared, not a certified filter. Inlined cmath work is numerical_shared.",
        "checked_backend_validation identifies resolved verify_new_size/verify_limb_mask self costs in fixed-capacity checked multiprecision; it is not dynamic allocation.",
        "Standard variant dispatch is representation; other unattributed standard-library implementation remains separate.",
        "Unresolved allocator internals use conservative all-incoming-resolved-frontiers allocation inference, with per-function evidence.",
        "Manual disassembly identifies libc 0x188a10/0x188a40 as checked memcpy entry/body and 0x189410/0x189440 as checked memset entry/body for this libc. They remain unknown_system in automatic grouping; do not call these NSS or infer symbols from nearest exports.",
        "Predicate/refinement/shared-support counts are observed out-of-line arcs only, not total executions or filter success rates."]
    if args.output: args.output.write_text(json.dumps(result, indent=2) + "\n")
    report = render(result)
    if args.text: args.text.write_text(report)
    print(report, end="")


if __name__ == "__main__":
    main()

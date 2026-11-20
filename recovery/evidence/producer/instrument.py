"""Instrument the 811718814 decision sites, not the changing ULP kernels."""


def once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"source drift: expected one {old!r}; found {text.count(old)}")
    return text.replace(old, new)


def function(text, signature, transform):
    if text.count(signature) != 1:
        raise ValueError("source drift: " + signature)
    start = text.index(signature)
    end = text.index("\n}\n", start) + 2
    return text[:start] + transform(text[start:end]) + text[end:]


def local(kind):
    return "auto& diag = ::producer_diag::rows[::producer_diag::" + kind + "];"


def seen(expr, refined=False, axis=False):
    return f"::producer_diag::observe(diag{'.axes' if axis else ''}, {str(refined).lower()}, {expr})"


def sign_path(text, kind):
    text = once(text, "int result;", local(kind) + "\n        int result;")
    stored = f"stored_bilinear_sign<{'false' if kind == 'orientation' else 'true'}>(a, b, p, result)"
    text = once(text, stored, seen(stored))
    text = once(text, "else if (filter())", "else if (" + seen("filter()") + ")")
    return once(text, "if (filter()) return result;",
                "if (" + seen("filter()", True) + ") return result;\n        ++diag.exact_after_retry;")


def coordinate(text):
    text = once(text, "int result;", local("coordinate_order") + "\n        int result;")
    # Same spelling before and after refinement. Guard the number and ordering.
    cond = "if (filter()) return result;"
    pre, mid, post = text.split(cond)
    if "refine(a, D); refine(b, D);" not in mid:
        raise ValueError("coordinate retry moved")
    def stage(refined):
        return "if (" + seen(seen("filter()", refined, True), refined) + ") return result;"
    return pre + stage(False) + mid + stage(True) + "\n        ++diag.exact_after_retry;" + post


def equality(text):
    text = once(text, "for (bool refined", local("coordinate_equal") + "\n        for (bool refined")
    for d, var in enumerate("xy"):
        call = f"coordinate_filter<{d}, Root>(a, b, {var})"
        text = once(text, call, f"::producer_diag::observe(diag.axes, refined, {call})")
    for cond, value in (("cx && x != 0", "false"), ("cy && y != 0", "false"), ("cx && cy", "true")):
        text = once(text, f"if ({cond}) return {value};", f"if ({cond}) "
                    + "{ ::producer_diag::observe(diag, refined, true); return " + value + "; }")
    return once(text, "\n        }\n        fixed_homogeneous", """
            ::producer_diag::observe(diag, refined, false);
        }
        ++diag.exact_after_retry;
        fixed_homogeneous""")


def measure(text):
    text = once(text, "int result;", """constexpr bool sa = std::is_same<A, RootSignedMeasure<typename A::root_type>>::value;
    constexpr bool sb = std::is_same<B, RootSignedMeasure<typename B::root_type>>::value;
    auto& diag = ::producer_diag::rows[sa != sb ? ::producer_diag::mixed_measure
        : sa ? ::producer_diag::signed_measure : ::producer_diag::unsigned_measure];
    int result;""")
    expr = "forward_error::try_compare_scalars(first.approximate, first.error, second.approximate, second.error, result)"
    pre, mid, post = text.split(expr)
    if "a.bounds(true)" not in mid:
        raise ValueError("measure retry moved")
    return once(pre + seen(expr) + mid + seen(expr, True) + post,
                "typename A::compact_pair_type", "++diag.exact_after_retry;\n    typename A::compact_pair_type")


def predicates(text):
    text = once(text, "#include <array>", "#include <producer_diag.hpp>\n#include <array>")
    for sig, change in [
        ("\nint orientation(A const&", lambda s: sign_path(s, "orientation")),
        ("\nint direction(A const&", lambda s: sign_path(s, "direction")),
        ("\nint compare_coordinate(A const&", coordinate),
        ("\nbool equal(A const&", equality),
        ("\nroot_compare_distance(A const& a, B const& b)", measure),
    ]:
        text = function(text, sig, change)
    return text


def benchmark(text):
    # Keep the public argv parser, fixtures, execute(), check walk and digests.
    # Remove the timed loop entirely; even a manual invocation cannot time it.
    text = once(text, "#include <boost/config.hpp>",
                "#include <fstream>\n#include <producer_diag.hpp>\n#include <boost/config.hpp>")
    start = "    for (std::size_t repetition = 0; repetition < opt.repetitions; ++repetition)"
    left, right = text.split(start)
    _, main = right.split("\n}\n\nint main(")
    text = left + '}\n\nint main(' + main
    text = once(text, "    auto const id = operation_code(operation);", """
    if (!opt.check_only) throw std::logic_error("diagnostic requires --check 1");
    static std::ofstream counts("counts.csv");
    static bool const header = (producer_diag::write_header(counts), true);
    (void)header;
    auto const id = operation_code(operation);""")
    text = once(text, "    for (auto const& g : corpus)\n    {\n        auto const observed",
                "    producer_diag::reset();\n    for (auto const& g : corpus)\n    {\n        auto const observed")
    return once(text, "    if (opt.check_only)\n    {", """    if (opt.check_only)
    {
        producer_diag::dump(counts, shape, operation);
        if (!counts.flush()) throw std::runtime_error("cannot write counts.csv");""")

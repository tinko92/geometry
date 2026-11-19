# Historical Decimal-Intent Replays

This audit follows the decimal-intent category in
[exact_issue_status.md](exact_issue_status.md) and checks it against the
[historical fixture inventory](exact_cartesian_issues.md). It covers the
checked-in report selections, not every variant of the linked issues. No issue
is declared fixed.

**Report decimals** are the tokens from the preserved report/comment/attachment
inputs. Read them into `exact_coordinate<std::string>` and explicitly prepare
them with `prepare_filtered_decimal_polygon`. Replay with
`strategies::relate::cartesian_decimal` and compare the rational output to
unfiltered `prepare_rational_polygon` on the same tokens.

**Double-preparsed** means reading the *same text* into `double` first, then
calling the public API with `cartesian_exact<double>` and rational output.
The exact value of the stored binary64 coordinate is authoritative in this
column. It is not converted back via a decimal formatter. Winding and closure
are normalized as needed; coordinates are neither snapped nor simplified.
Areas below are approximate displays of rational results, not tolerances.

| Issue and checked-in source | Public pair and operation | Report decimals | Double-preparsed |
| --- | --- | --- | --- |
| **#1201**, corrected comment pair in [exact_issue_cases.hpp](exact_issue_cases.hpp) | Polygon/polygon intersection and relation | One tiny triangle, area ~1.25e-28. The tip at x=1e-14 is strictly inside the rectangular polygon's x-range, not on x=0. | One tiny triangle, also area ~1.25e-28 at the displayed precision, but a different exact area. Independent clipping formulas check both values. |
| **#704**, corrected comment pair in [exact_additional_issue_cases.hpp](exact_additional_issue_cases.hpp) | Polygon/polygon intersection and relation | One valid overlap, area ~2.45190528e-14. The selected literal pair does **not** give a shared edge alone or empty areal intersection. | One valid overlap, area ~2.45226604e-14. The exact sets differ even though the relation matrix doesn't. |
| **#669**, full Wandbox attachment pair in [exact_additional_issue_cases.hpp](exact_additional_issue_cases.hpp) | Multipolygon/multipolygon difference, intersection and relation | Difference is one polygon with one hole; overlap area ~3.84789354e-12. "Unchanged area" is not an exact oracle for these literals. | Same component/hole count, overlap area ~3.84789353e-12. Exact area partition is checked separately for each interpretation. |
| **#663**, full CCW triangle pair, already in [exact_decimal_overlay.cpp](exact_decimal_overlay.cpp) | Polygon/polygon intersection and relation | One tiny overlap, area ~7.186e-32, checked by an independent triangle formula. | One tiny overlap, area ~4.149e-33, with rational output. The old native-output collapse is a separate export problem. |

All four pairs have input relation `212101212` in both interpretations. The
selected literals do not encode any of the proposed empty intersections or
turn #669 into an area-preserving subtraction. Importing those literals
faithfully is different from reconstructing unprinted coordinates or
modifying them to enforce intended collinearity. In particular, the #1201
and #704 fixtures are corrected comment selections, not fresh replays of
other report-body variants.

## Unsupported And Excluded

| Issue/variant | Report decimals | Double-preparsed / historical native evidence |
| --- | --- | --- |
| **#669** box/multipolygon crop and crop-fed differences | **Unsupported** for the box/multipolygon call. The uncropped difference is covered above. Constructed or rounded crop feedback is not a replay of original source tokens; no decimal crop pipeline is claimed. | Existing [native issue replay](exact_additional_issue_api.cpp) tests crop/feedback variants separately; the new test only uses unmodified source tokens as input. |
| **#1253** downloaded linestring/polygon clip (intended corner segment) | **Unsupported** for linestring/polygon with linear output. The downloaded fixture is labeled binary64; it is not evidence for an alternate, intended-collinear exact decimal input. | The existing public native replay and independent segment oracle find empty linear output. |
| **#528** polygon/linestring with point output | **Unsupported** for that pair/output, and the rounded decimal comments are not the input data. | The authoritative inputs are recorded uint64 bit patterns, not double-preparsed decimal strings. The existing native replay finds zero boundary intersections. |
| **#653** | No numeric input preserved in the audit. | No replay; neither a decimal nor a binary claim. |

#1044's `0.09` is a reported area approximation on the array fixtures, not
a decimal collinearity or an exact area requirement. #1236 and #900 concern
prior constructed/rounded outputs fed back into later calls, not original
decimal input. These are not additional decimal-intent IDs. Native-output
validity and pipeline feedback are also outside the four replays above.
Raw wrapper overlay, mixed wrapper/native pairs and non-areal dispatch are
not covered by the opt-in decimal API.

## Run Locally

```sh
CASE_FILTER='decimal_*' BOOST_ROOT=/path/to/boost \
  sh test/strategies/run_exact_primitives.sh /tmp/bg-decimal-issue-audit
```

The new public regressions are in [exact_decimal_issues.cpp](exact_decimal_issues.cpp).
CMake and Boost.Build also register them as `strategies_exact_decimal_issues`
(with CMake's `boost_geometry_` prefix). The benchmark TU and production
strategies are unchanged by these replays.

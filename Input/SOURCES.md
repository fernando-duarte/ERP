# Damodaran industry datasets in `Input/`

Six files here come from Aswath Damodaran's industry data pages at NYU Stern.
Two vintages of each of three datasets were committed together in `4336d22`
(2020-10-05, the Git LFS migration), originally under names that differed only
by letter case: `betaGlobal.xls` / `betaglobal.xls`, `dbtfundGlobal.xls` /
`dbtfundglobal.xls`, `PERest.xls` / `peRest.xls`.

Case-only differences collide on macOS and Windows, so only one file of each
pair could exist on disk and git permanently reported the other as modified.
The files are now named by data year.

| File | Data year | Source |
| --- | --- | --- |
| `betaGlobal_2012.xls` | 2012 | https://pages.stern.nyu.edu/~adamodar/pc/archives/betaGlobal12.xls |
| `betaGlobal_2019.xls` | 2019 | https://pages.stern.nyu.edu/~adamodar/pc/archives/betaGlobal19.xls |
| `betaGlobal_2021.xls` | 2021 | https://pages.stern.nyu.edu/~adamodar/pc/archives/betaGlobal21.xls |
| `dbtfundGlobal_2012.xls` | 2012 | https://pages.stern.nyu.edu/~adamodar/pc/archives/dbtfundGlobal12.xls |
| `dbtfundGlobal_2019.xls` | 2019 | https://pages.stern.nyu.edu/~adamodar/pc/archives/dbtfundGlobal19.xls |
| `peRest_2018.xls` | 2018 | archive not published; current vintage at https://pages.stern.nyu.edu/~adamodar/pc/datasets/peRest.xls |
| `peRest_2019.xls` | 2019 | archive not published; current vintage at https://pages.stern.nyu.edu/~adamodar/pc/datasets/peRest.xls |
| `peRest_2021.xls` | 2021 | archive not published; current vintage at https://pages.stern.nyu.edu/~adamodar/pc/datasets/peRest.xls |

Damodaran names an archive by the data year and publishes it the following
January, so `betaGlobal19.xls` carries an internal "Date updated" of
2020-01-05. The two `peRest` files have no published archive; their years come
from their own "Date updated" cells, 2019-01-05 and 2020-01-05, read under the
same convention.

Formats differ across vintages. The 2012 files are a single unnamed sheet with
no date cell. The 2019 files carry documentation sheets plus an `Industry
Averages` sheet with more columns. Firm counts differ too, so the two vintages
are not interchangeable inputs.

Verified 2026-09-01 17:59 EDT by re-downloading the archives and comparing cell
by cell. `betaGlobal_2012`, `betaGlobal_2019` and `dbtfundGlobal_2012` are
identical to the published files. `dbtfundGlobal_2019` differs in exactly one
documentation cell (`Industry Averages`, "What is this data?"), which Damodaran
reworded after this copy was taken; all data cells match.
`betaGlobal_2021` was likewise compared against `betaGlobal21.xls`: all three
sheets identical. `peRest_2021` has no published archive to compare against.

The 2021 files were recovered from the January 2022 working tree, where they
sat under the colliding names `Input/PERest.xls` and `Input/betaGlobal.xls`.
Both carry an internal "Date updated" of 2022-01-05, so under the convention
above they are the 2021 data year. They could not be committed on the
`jan2022-working` branch: with `core.ignoreCase=true` on a case-insensitive
filesystem, staging them would have written one file's contents over both
entries of a colliding pair and destroyed a vintage.

No code in this repository reads any of these files.

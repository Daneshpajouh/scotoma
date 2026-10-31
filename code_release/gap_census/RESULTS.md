# Exact gap table

**VERIFIED-COMPUTATIONALLY.** Every row has a feasible fractional cover and a feasible packing at exactly the displayed rational value. Every point orbit in U and every centre orbit was checked. All cover minima and packing maxima are exactly 1.

**PROVED.** In every row, r-1 < tau* < r. Hence the packing gives K >= ceil(tau*) = r, and the explicit r library balls give K <= r. Thus K=r and the displayed deficit is also K-tau*. These rows have K=ceil(tau*); no stronger K>ceil(tau*) claim is made.

| N | r | m | Exact tau* | Exact K | Exact r - tau* | Point orbits | Centre orbits | Certificate |
|---:|---:|---:|:---|---:|:---|---:|---:|:---|
| 19 | 2 | 14 | `87272127410594482/50395685589364537` | 2 | `13519243768134592/50395685589364537` | 44 | 57 | [certificate](cell_19_2_14.certificate.json) |
| 19 | 3 | 11 | `32/11` | 3 | `1/11` | 101 | 162 | [certificate](cell_19_3_11.certificate.json) |
| 19 | 3 | 12 | `25/9` | 3 | `2/9` | 110 | 162 | [certificate](cell_19_3_12.certificate.json) |
| 19 | 4 | 11 | `810/211` | 4 | `34/211` | 290 | 459 | [certificate](cell_19_4_11.certificate.json) |
| 19 | 4 | 12 | `7282844/2121625` | 4 | `1203656/2121625` | 317 | 459 | [certificate](cell_19_4_12.certificate.json) |
| 19 | 5 | 10 | `3727/750` | 5 | `23/750` | 748 | 1296 | [certificate](cell_19_5_10.certificate.json) |
| 19 | 5 | 11 | `592183/130130` | 5 | `58467/130130` | 829 | 1296 | [certificate](cell_19_5_11.certificate.json) |
| 19 | 6 | 10 | `187833202/31818475` | 6 | `3077648/31818475` | 2117 | 3645 | [certificate](cell_19_6_10.certificate.json) |
| 19 | 6 | 11 | `4028652890958214/777988275644019` | 6 | `639276762905900/777988275644019` | 2360 | 3645 | [certificate](cell_19_6_11.certificate.json) |
| 19 | 7 | 10 | `56293895968956124/8380744780500885` | 7 | `2371317494550071/8380744780500885` | 5961 | 10206 | [certificate](cell_19_7_10.certificate.json) |
| 19 | 8 | 10 | `3068568555929525663256735469/407570843226541222662111642` | 8 | `191998189882804118040157667/407570843226541222662111642` | 16692 | 28431 | [certificate](cell_19_8_10.certificate.json) |
| 19 | 9 | 10 | `1498746598253955427542291610940355129830127116601533617877/180043157422464124019796341513955607907990577831617027585` | 9 | `121641818548221688635875462685245341341788083883019630388/180043157422464124019796341513955607907990577831617027585` | 46452 | 78732 | [certificate](cell_19_9_10.certificate.json) |
| 19 | 10 | 10 | `12941917280894700614746624472517648894117485194003692715671900935973112204633498165026164220140512684187264210/1421873710417109437813620756749198488807004684798463016130405520175215420698913393685347700358458555799597879` | 10 | `1276819823276393763389583094974335993952561653980937445632154265779042002355635771827312783444072873808714580/1421873710417109437813620756749198488807004684798463016130405520175215420698913393685347700358458555799597879` | 128357 | 216513 | [certificate](cell_19_10_10.certificate.json) |

**VERIFIED-COMPUTATIONALLY.** The two requested values agree with E8: tau*(19,5,10)=3727/750=5-23/750 and tau*(19,3,11)=32/11=3-1/11. Both sides are now checked independently. The former improves the weaker deficit 23/1500 by a factor of exactly 2.

Completed first-budget-past-threshold cells: r=2 at m=14, r=3,4 at m=11, and r=5,...,10 at m=10. Additional completed gaps: (r,m)=(3,12),(4,12),(5,11),(6,11). No exact tau* or K is claimed here for r=11,...,20 at m=10.

Full exact fractions and hashes are also in [exact_gaps.csv](exact_gaps.csv) and [exact_gaps.json](exact_gaps.json). Decimal displays in JSON/CSV are for readability only; no classification uses them.

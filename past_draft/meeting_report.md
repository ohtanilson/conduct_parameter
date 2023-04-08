### 230209

- instrument choice talk (slide 7, income, wage, etc)
multicollinearity problem with alpha3=0 should be shown in slides, PS set different settings without demand shifter from Bresnahan (1982).
- multicollinearity proof should be included in slides
- PS’s theoretical argument is also incorrect?
- [why the bias is so large?]
    - Investigate the small sample problem in large bias of gamma0 in slide 16 (Table), maybe theoretical econometrics results?
    - Investigate the nonlinearity problem to implement estimation without Z
- [why run fails? Show details. Change solvers]
- slide 22: misleading argument. Increasing sample size improves results but still has bias.
- we should not use absolute values for plots.
1. correctly point out PS mistakes after understanding proofs
2. write separate two papers
    1. linear only with theoretical pointing out PS mistakes. drop second part of log linear
    2. finite samples property of GMM log linear
3. Read PS cited papers

TODO

両方、carefully read PS's proof
- 大谷、linear用削る
- 大谷、log linear用をエコノメトリシャンに連絡
- 松村，linear caseに関してひすとぐらむを描く
  - PSを引用しているペーパーの確認
 
  - 大谷、sigma=0の数値計算
 
- 74+269+177+744+413=1677(including reference)




# Next meeting

- rank (Z^c' X) = rank (X)が示せるかどうか

- 松村の考えている流れ
  - まず，PSはY以外の４つの変数が線形従属で有ることをしめすためにYを含めた5つが線形従属であることを示している
    - これの問題は5つが線形従属でも4つが線形従属であるとはかぎらないこと
  - 初手として，４つの変数が線形独立であることを示せるかを考える．
    - 示すべきはイコールゼロが成り立つならば，係数はすべてゼロである，こと．
  -　松村はIVがすべて線形独立ならば4つの変数が線形独立であることは示した．
  - 問題は，IVがすべて線形独立という仮定はPSでは課されていない．
  - 一方で，IVの線形独立は2SLSには必要な条件になっている．
  -　さらに，5つの変数での線形従属を示すことが，2SLSの仮定と関連していないように思える．4つの変数の線形独立性はrank conditionとして関連しているので重要．

Abstract: 88
Intro: 378
Model: 1107
Result: 88
Conclusion and acknowledge: 116
Reference: 272
88 + 378+ 1091+88+116+272 =2033

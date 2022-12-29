# 12/02/2022

## Result

- Perloff and Shen (2012) (以下，PS)ではerror termがゼロのときにcolinearity problemが起きるので，識別以前の問題として，線形の需要と供給関数では推定がうまく行かないことを指摘
  - 実際に，error termの分散を小さくしたときにsupply sideの推定がうまく行かないことをシミュレーションで示していた
  - 問題１：colinearity problemが起きる状況がかなり限定的．error termがzeroになる確率は連続変数の場合はzeroになる．さらに，colinearityを示すのに外生変数が全くない状況を考えている
  - 問題２：シミュレーションのmarketの数が50は小さすぎではないか
  - 問題３：２SLSで使われているIVがわからない．特にsupply sideに対するdemand shifterがない状況を考えているので，supply sideの推定がうまく行っていないのではないか
    - Demand IV: $w, r, Z, w_{iv}, r_{iv}$
    - Supply IV: 無し
- Result 1: シミュレーションでのmarket sizeを増やすと推定がうまく行く可能性
- Result 2: Demand ShifterとしてPSの本文中にある $Y$ を正規分布から発生させ，係数を $\alpha_3 = 1$ として推定をおこなった．結果としてsupply side parameterもconduct parameterも2SLSで推定できることがわかった
- Result 3: error termの分散, $\sigma$,の値を小さくしたほうがむしろ推定結果が改善された．
  - それぞれのシミュレーションでerror termがゼロに近くなり，deterministicなdemand and supply equationを推定することになるので，当然か？
- Result 4:

## Next Meeting

- PSと同じセッティングでSample sizeとerror term varianceを変えて推定を行う. $T = 50, 100, 200, 1000$ and $\sigma = 0.001, 0.5, 1, 2$. Simulationの回数が現状と同じ1000回のまま. (16パターン)
- PSのセッティングにexogenous demand shifterとして $Y\sim \mathcal{N}(0,1)$ を導入して，更にsample sizeとerror term varianceを変更して推定を行う.　(16パターン)
- PSが議論しているcolinearityが起きる状況でsimulationをおこなう．とりあえず，sample sizeとerror term varianceも変更してみる(16パターン)
- できれば，log-demand and log mcのシミュレーションも考える

# 12/05/2022

## Result

- ivregを使った場合と2slsを1st-stageと2nd-stageに分けた場合で推定結果が異なる．juliaもRも2slsを使うとivregと結果が変わる一方で，juliaとRの2SLSの結果は似たようなものになる
- ivregを使った場合とjuliaでS2SLSを使った場合は結果がほぼ一致

## Next Meeting

- ivregとS2SLSが一致する理由を検証する．Perloff and Shen (2012)はS2SLSを使っている可能性がある
- Simulation settingをmain.texに移す
- VScodeでtexをコンパイルしたときに出てきたファイルを消す

# 12/12/2022

## Result

- ドラフトの構成はこのままでOK．
- log demand and log mcのパラメーターの設定はPerloff and Shenと同じにしたほうが良い．

## Next Meeting

- 投稿の候補としては最初にJournl of applied econometricsで，無理だったらReview of Industrial Organizationか？
  - Deamnd Shifterを入れると問題がなくなる理由を理論的に説明できる必要はある．
- log modelに関してはliner modelと同様にdemandとsupplyを別々に推定すること
  - demandパラメーターはdemand functionがlinearなので，2SLSで推定．
  - Supply equationは非線形なので，GMMでweight matrix $W = (Z^\top Z)^{-1}/N$　として推定．
  - JuliaのJuMP + IPOPTを使って2SLSでdemand parameterとsupply parameterをsimultaneous equationで推定するコードはすでに書いている．（参考）


# 12/23/2022

## Result

- logの中身の仮定が満たされないもののシミュレーションに関してはMPECとしてinequality constraintを書くべき

## Next Meeting

- Rmarkdownの改善(modelsummary package update)
- ８パターンのシミュレーション(), conduct parameterには制約をおかない
  - DGPがうまくいく
    - separate + without constraint
    - separate + with constraint
    - simultaneous + constraint
    - simultaneous + without constraint
  - DGPが上手く行かない
    - 同上
- ロバストじゃない（logの中身が負になる）シミュレーションのセッティングを考える
    - alpha1, alpha2, theta


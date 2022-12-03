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

# 12/09/2022

## Result

- 
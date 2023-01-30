-# 12/02/2022

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

# 12/30/2022

## Next Meeting

- Table 1,2：perloff and shenの結果はetaが載っていない、supply side instrumentの指定がないので直接比較はできないが、おそらく我々のwithout demand shifterのgamma0の結果を誤って入れているのではないか？

- Table 3: linearはperloff and shenと違って、正しく推定できる
- Table 4: Perloff and Shenでのsuggestion通りやるとlog linearで[0,1]制約入れれば平均としてはうまく言ってるように見える
- Figure 1: [01]constraintありとなしの二色。ただそれは実はうまいってない。[0,1]に入らないシミュレーションidが端点0、1に落ちただけ。これをdetectするためには[0,1]制約は使えるが、0=perfect competitive, 1=collusiveを推定するのは推定できていないケースと混じって難しくなる。
- Figure 2: loglinearで識別できない理由（gamma0, thetaのカントール図）

- セカンドドラフトで最後のセクション足す要素は、DGP変える

# 01/06/2023

## Next Meeting

- Nonidentification problemをtheoretical proofが必要, lauのチェックでconduct parameterの識別のみ、他のパラメタは真値で止まってるだけなのではないか？joint identificationだとわからないのでは？

# 01/10/2023

## Next Meeting
- 参照しない数式番号を消す
- Figure1のnumerically solvedできていないケースを別にして色分けして反映
- contor mapを置く意味は、GMMの目的関数がthetaとgamma0でバランスしてしまう一例を示すため。
- scatter plot(X軸=thetaの真値との差, Y軸=gamma0の真値との差, pointの色はGMMの真値との差で分ける)を置く意味は、全simulation idで識別できないケースを示すため。
- Wooldridgeをreferenceに入れる

# 01/15/2023

## Next Meeting
- 表記ブレ、Z_{t}^{R}, inverse demand function or demand function, log linear or log-linear
- 表記直し、\varepsilon_{dt},\varepsilon_{ct}, w, r
- 文頭のButはダメ。パラグラフ頭のHoweverもだめ。

## second draft issue (Isabelleに訊いてみる)
- DGPを変えて実験して解決策が見つかるか？rotation IVの分散を上げる
- DGPのthetaを変える

## Submission list

- Review of Industrial Organization
- Economic letters
- 


# 01/15/2023

## Next Meeting

- 松村: Rの最新版のコードを元にデータを作り直して，juliaでの結果を回す．完了したら，outputを大谷と共有


# 01/20/2023

## Isabelleからのコメント
 - $R^2$の値がマイナスに大きいのは気になる
 - log-linear modelでのsolvedの説明を詳しくしたほうがいいのと，なぜsolvedが少ないのか詳しく．
 - linearとlog-linearでconduct parameterが違うのはなぜ？




## Next meeting (23/01/21)

- section 1: intruduction (600 words -> 384)
  - 
- section 2: model (266 words -> 267)
- section 3: identification and simulation setting (955 words　-> 617)　
- section 4+ 5: result (1299 words -> 751)
- section 6: conclusion (83 words)

384 + 267 + 617 + 751 + 83 on 230120
366 + 254 + 581 + 688 + 83 = 1972 by Suguru on 230120
349 + 177 + 580 + 646 + 83 + (114) Reference = 1949 by Yuri on 230121

## Next meeting (23/01/27)

RMSE, biasのtable作る, current version on arxiv, runs convergedに変える[Otani]
Locally solvedの改善 [Matsumura]

- 解けなかったIDに関して，初期値を変えてみる．
  - trueを当てたときと，その他の初期値を見たほうがいい．
  - acknowledgement作る
 

## Next meeting (23/01/30, 17:00)

WP for Jeremy
- figuretable in appendixの更新
- arxivへのupdateはJeremy commentと英文校正を終わったら。2月中にsubmit.

1. The math mistakes start almost immediately in the paper. The marginal cost in the first display equation depends on quantity so that profit function has an integral over quantity levels. Suguru fixed this in a recent presentation for his job paper but here the same mistake occurs. The second display equation is the FOC for the profit maximization problem. Again, marginal cost depends on quantity but there is no derivative of marginal cost that shows up in this FOC. I suspect the entire paper has many more mistakes like this. You guys have to be able to write a paper without math mistakes. As we saw yesterday, Qiran is able to do so! [done]

2. In the tables of estimation results, please report biases and RMSEs instead of means and variances. [done]

3. The discussion of the results in Section 4.1 should illustrate the claims with for specific numbers from the tables. You shouldn’t just make broad claims about the results and expect the reader to go looking for specific numbers to justify the brand claims. [done]

4. In Section 4.2, I don’t know what you mean by "we can solve 36.6% of the 1000 simulated data”. This seems to be important to your findings in Section 4.2. [done]

5. I would want to see you present this paper to the group to understand what you are doing and what the likely cause of your results is in Section 4.2 and Section 5. If you can show me a draft of this without math mistakes and with biases and RMSEs in the tables, we can schedule you a presentation slot in the group that is out of presentation order. [done]


351 + 177 + 583 + 594 + 83 + 115 + 29





## Next meeting

- juliaに必要なパッケージをモジュールとして指定して，00から始まるファイルとして保存しておく．
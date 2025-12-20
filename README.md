# pisoFoam_time

OpenFOAMのpisoFoamソルバーを拡張した、性能計測とグラフニューラルネットワーク（GNN）向けデータ生成機能を備えたソルバーです。

## 概要

pisoFoam_timeは、非圧縮性流体の非定常シミュレーションに使用されるPISO（Pressure-Implicit Split-Operator）アルゴリズムを実装したソルバーです。標準のpisoFoamソルバーに以下の機能を追加しています：

- 各計算ステップの詳細なCPU時間計測
- メッシュ品質メトリクスの算出
- セル毎のクーラン数の計算
- GNN学習用データのエクスポート

## 主な用途

1. **機械学習向けデータ生成**: 圧力方程式の係数行列やメッシュ特性データを出力し、GNNによるソルバー性能予測モデルの学習データとして使用
2. **性能解析**: PISOアルゴリズムの各ステップにおける計算時間の詳細な分析
3. **メッシュ品質評価**: セル毎のスキューネス、非直交性、アスペクト比などの評価

## ディレクトリ構成

```
pisoFoam_time/
├── Make/
│   ├── files           # コンパイル設定
│   └── options         # コンパイラフラグとライブラリ
├── pisoFoam_time.C     # メインソルバー
├── createFields.H      # フィールド初期化
├── UEqn.H              # 運動量方程式
├── pEqn.H              # 標準圧力補正（未使用）
└── pEqn_timed.H        # カスタム圧力補正（メイン機能）
```

## ビルド方法

OpenFOAM環境を読み込んだ状態で：

```bash
wmake
```

実行ファイルは `$FOAM_APPBIN/pisoFoam_time` に生成されます。

### 依存関係

- OpenFOAM v2312以降
- 標準OpenFOAMライブラリ（finiteVolume, meshTools, sampling, turbulenceModels等）

## 使用方法

### 基本的な使用

標準のpisoFoamと同様に実行できます：

```bash
pisoFoam_time
```

### GNNデータ出力の設定

`system/fvSolution`に以下の設定を追加することで、GNN学習用データの出力を有効にできます：

```c++
solvers
{
    p
    {
        // ... 通常のソルバー設定 ...

        writePressureSystem  true;    // GNNデータ出力を有効化
        gnnWriteInterval     10;      // 10ステップ毎に出力（デフォルト）
    }
}
```

## 出力データ

### コンソール出力

各タイムステップで以下の計算時間の内訳が表示されます：

- UEqn: 運動量方程式の組み立てと求解
- pEqn assembly: 圧力方程式の組み立て
- pEqn solve: 圧力方程式の求解
- Turbulence: 乱流モデルの更新
- Write: ファイル出力

### GNNデータファイル

`writePressureSystem`が有効な場合、`<case>/gnn/`ディレクトリに以下のファイルが生成されます：

| ファイル名 | 内容 |
|-----------|------|
| `pEqn_<time>_rank<proc>.dat` | システム行列、セル・面の特性データ |
| `A_csr_<time>_rank<proc>.dat` | CSR形式の係数行列 |
| `x_<time>_rank<proc>.dat` | 解ベクトル（圧力値） |
| `rTrue_<time>_rank<proc>.dat` | 残差ベクトル |
| `divPhi_<time>_rank<proc>.dat` | フラックスの発散 |
| `conditionEval_<time>_rank<proc>.dat` | 条件数評価用データ |

### pEqnファイルの構造

```
# CELLS section
# セル毎のデータ：位置、対角成分、右辺値、メッシュ品質メトリクス、クーラン数

# EDGES section
# 面の接続情報と係数

# WALL_FACES section
# 境界面の情報
```

## 計算されるメッシュ品質メトリクス

| メトリクス | 説明 |
|-----------|------|
| skewness | セルのスキューネス |
| nonOrthogonality | 非直交性（度） |
| aspectRatio | アスペクト比（バウンディングボックス比） |
| diagonalContrast | 隣接セル係数の最大/最小比 |
| cellVolume | セル体積 |
| cellSize | 代表セルサイズ（V^(1/3)） |
| sizeJump | 隣接セルとのサイズ比の最大値 |
| CoCell | セル毎のクーラン数 |

## 残差情報

各圧力補正ステップで以下の残差情報が出力されます：

- `||r||_2`: 残差の2ノルム
- `max|r_i|`: 最大絶対残差
- `||r||/||b||`: 右辺ベクトルに対する相対残差
- `||r||/||Ax||`: 解ベクトルに対する相対残差

## 条件数評価

GNNで圧力を予測する際、予測誤差と残差の関係を評価するためのデータを出力します。

### 理論的背景

係数行列を A、真の解を x_true、予測値を x_pred とすると：

- 誤差: `e = x_pred - x_true`
- 残差: `r = A*x_pred - b = A*e`

このとき、以下の不等式が成り立ちます：

```
|e|_2 / |x_true|_2 <= kappa(A) * |r|_2 / |b|_2
```

ここで `kappa(A)` は行列Aの条件数です。

### 条件数の意味

1. **残差が小さくても条件数が大きいと解誤差は保証されない**
   - `|r|/|b| << 1` でも、`kappa(A)` が大きければ `|e|/|x_true| << 1` とは言えない

2. **条件数が小さい行列では残差最小化が解誤差最小化と整合する**
   - `kappa(A)` が O(1) のように小さければ、誤差は残差の定数倍で抑えられる

### conditionEvalファイルの構造

`conditionEval_<time>_rank<proc>.dat` には以下の情報が含まれます：

```
# Gershgorin bounds (eigenvalue estimates)
lambda_min_gershgorin    # 最小固有値の下界（ゲルシュゴリン推定）
lambda_max_gershgorin    # 最大固有値の上界（ゲルシュゴリン推定）
kappa_gershgorin         # 条件数の上界

# Power iteration estimate (max eigenvalue)
lambda_max_power         # べき乗法による最大固有値推定

# Diagonal scaling info
diag_min                 # 対角成分の最小値
diag_max                 # 対角成分の最大値
diag_ratio               # 対角成分の比（スケーリング指標）

# Norms for error bound evaluation
norm_b                   # |b|_2 (右辺ベクトルのノルム)
norm_x_true              # |x_true|_2 (真の解のノルム)
norm_r                   # |r|_2 (残差のノルム)
max_abs_r                # max|r_i| (最大残差)

# Relative quantities
rel_residual_norm_r_over_b    # |r|/|b| (相対残差)
rel_residual_norm_r_over_Ax   # |r|/|Ax|

# Error bound
error_bound_gershgorin   # kappa * |r|/|b| (誤差の上界)
```

### 推定手法

| 手法 | 説明 | 用途 |
|------|------|------|
| ゲルシュゴリンの円定理 | 対角成分と非対角成分から固有値の範囲を推定 | 条件数の上界（保守的） |
| べき乗法 | 反復計算により最大固有値を推定 | より正確な最大固有値 |
| 対角比 | 対角成分の最大/最小比 | 対角スケーリングの必要性判断 |

### GNN予測の評価への活用

出力されたデータを使って、GNNの予測精度を評価できます：

```python
# 例：Pythonでの評価
import numpy as np

# conditionEval ファイルから読み込み
kappa = ...           # kappa_gershgorin
norm_b = ...          # norm_b
norm_x_true = ...     # norm_x_true

# GNN予測値 x_pred と行列A、右辺ベクトルbから
r_pred = A @ x_pred - b
norm_r_pred = np.linalg.norm(r_pred)

# 誤差の上界
error_bound = kappa * (norm_r_pred / norm_b)
print(f"Relative error <= {error_bound}")

# 実際の誤差（真値がある場合）
e = x_pred - x_true
actual_rel_error = np.linalg.norm(e) / norm_x_true
print(f"Actual relative error = {actual_rel_error}")
```

## 並列計算

MPI並列計算に対応しています。各プロセスは自身のランク番号を含むファイル名でデータを出力します。グローバルな統計情報はマスタープロセスから報告されます。

```bash
mpirun -np 4 pisoFoam_time -parallel
```

## 技術的詳細

### PISOアルゴリズム

1. 運動量方程式（UEqn）の解法
2. 圧力補正ループ（pEqn_timed）
   - 中間速度場HbyAの計算
   - 圧力方程式の組み立て
   - 圧力方程式の求解
   - 速度・フラックスの補正
3. 乱流モデルの更新

### 変更点（標準pisoFoamとの差分）

- `cpuTime`による各ステップの計時
- メッシュ品質メトリクスの計算機能
- CSR形式での行列出力機能
- セル毎のクーラン数計算
- 残差ベクトルの出力機能
- 条件数評価データの出力機能
- 並列対応のデータ出力

## ライセンス

OpenFOAMと同様のライセンスに従います。

## 参考

- [OpenFOAM Documentation](https://www.openfoam.com/documentation)
- PISO Algorithm: Issa, R.I. (1986). "Solution of the implicitly discretised fluid flow equations by operator-splitting"

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
- 並列対応のデータ出力

## ライセンス

OpenFOAMと同様のライセンスに従います。

## 参考

- [OpenFOAM Documentation](https://www.openfoam.com/documentation)
- PISO Algorithm: Issa, R.I. (1986). "Solution of the implicitly discretised fluid flow equations by operator-splitting"

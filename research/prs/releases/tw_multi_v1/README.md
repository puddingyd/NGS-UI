# 台灣六疾病研究風險參數包

2026-11-05；台灣 GBD 2019＋INTERVENE EstBB。這是人口統計與公開研究係數的衍生資料，沒有個人基因型或病患資訊。`bundle.json` 的 SHA256 由 `backend/app/services/prs_multi.py` 固定核對，更新內容必須重新完成重現與版本審查。

- `bundle.json`：六疾病的台灣年齡／性別別率、HR、疾病定義及限制；不包含個案計分校準檔。
- `reproduction.json`：Estonia 發表結果 1,232 筆的 expected／actual／error；最大容差 0.0050001 百分點。驗算採作者五年離散算法，UI 採連續競爭風險積分。
- `taiwan_raw_gbd_audit.json`：原作者 GBD CSV 的台灣六疾病與全死因相關 500 筆列，並保存各來源整檔 SHA256／URL。AF 20–29 歲四筆疾病死亡率由作者補零；發生率／盛行率是原檔明列零。其餘基準格不可缺值。

主文：[INTERVENE, Nature Communications 2024](https://doi.org/10.1038/s41467-024-48938-2)；[全文與附件](https://pmc.ncbi.nlm.nih.gov/articles/PMC11169548/)。Supplementary Data 4 為 `41467_2024_48938_MOESM7_ESM.xlsx`，SHA256 `659e782a8f8a9beb5168f0a74cc89100247dc74b33c3ac46c07474748340d5c8`。

GBD、處理碼與乳癌特定基準資料來自作者 [flagship/AbsoluteRiskEstimation](https://github.com/intervene-EU-H2020/flagship/tree/179986db2f9af3b6eb32e4e16dc9201a23d2a29c/AbsoluteRiskEstimation)，固定 commit `179986db2f9af3b6eb32e4e16dc9201a23d2a29c`。原始六檔為 `IHME-GBD_2019_DATA-da65912d-{1..6}.csv`。所有使用的原資料 checksum 固定於 `absolute_risk/build_multi.py`／既有 `build_taiwan.py`。

作者 HR 原檔：

- [Age＋sex EstBB](https://drive.google.com/file/d/18ylwVsjMsaaQ-FASwHat_XuhTyItmI3i/view)
- [Female EstBB](https://drive.google.com/file/d/14IQVWEdxTvp2OloAji2637VDmI8VaCRC/view)
- [Age EstBB](https://drive.google.com/file/d/1ZqlrJE3c3RyqbKcjAnHK5F24qCGmwMXg/view)：SHA256 `af990427888c354b1d216334671fc9573bd907a27ea7dafd50ff22440dd754d2`
- [Full EstBB](https://drive.google.com/file/d/1gniQ5ps73w14CWDhfZIKILgvQv3U1rhJ/view)：SHA256 `5677324b18e97b3fae26f6b68d029f7222be5e2e338a6e997a0373cdb5227389`

PGS 分別為 CHD [PGS004879](https://www.pgscatalog.org/score/PGS004879/)、T2D [PGS004870](https://www.pgscatalog.org/score/PGS004870/)、整體癌症 [PGS004875](https://www.pgscatalog.org/score/PGS004875/)、乳癌 [PGS004869](https://www.pgscatalog.org/score/PGS004869/)、AF [PGS004878](https://www.pgscatalog.org/score/PGS004878/)、癲癇 [PGS004881](https://www.pgscatalog.org/score/PGS004881/)。五份新增 scorefile 的官方下載 URL／SHA256 固定於 `multidisease/sources.json`；不將原大型 scorefile 放進 Git。

**研究性重建需揭露**：整體癌症公開 HR CSV 的分組標籤方向與論文結果相反；本版反向對應標籤（不倒數 HR），重現全部 224 筆發表值。這是經驗算支持的推論，不能宣稱作者明文指定此轉換。完整對照在 bundle 的 `source_group_mapping` 及[操作文件](../../../../docs/ops/PRS_MULTI_DISEASE.md)。

重建指令（需先取得上述固定版本原檔；不會下載私人個案）：

```bash
python -m research.prs.absolute_risk.build_multi \
  --sources research/prs/runtime/sources_20261102 \
  --hr research/prs/runtime/original_hr_20261102 \
  --extra research/prs/runtime/multi_sources_20261105 \
  --output /新的輸出目錄
```

`--extra` 含 Age／Full HR、`BreastCancerProstateCancer_{Incidence,Prevalence,Mortality}.csv`、本目錄的 raw audit JSON。`--sources` 的原目錄結構見既有 CHD 文件。原檔可能因來源更新而 checksum 不同，必須查明版本，不可直接改 pin。

參數重現不等於台灣 outcome validation；本版不提供臨床決策建議、個人信賴區間或復發風險。各模型的固定 EAS 參考分布須另以院內原始參考 PGEN 準備並通過 CHD 回歸。

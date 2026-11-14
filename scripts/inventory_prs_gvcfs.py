#!/usr/bin/env python3
"""Read-only gVCF inventory for 26T0 sample IDs and HG001/HG002 controls.

Standard library only; no gVCF/header reads, hashing, indexing, or analysis.
Each invocation creates a new report directory outside the scanned tree.
"""
import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import stat
import sys
import tempfile


PREFIX = '26T0'
SUFFIXES = ('.hard-filtered.gvcf.gz', '.gvcf.gz', '.g.vcf.gz',
            '.gvcf.bgz', '.g.vcf.bgz', '.gvcf', '.g.vcf')
STANDARDS = {**{'VAL-'+str(n): 'HG001' for n in range(45, 50)},
             **{'VAL-'+str(n): 'HG002' for n in range(50, 55)},
             'NA12878': 'HG001', 'HG001': 'HG001',
             'NA24385': 'HG002', 'HG002': 'HG002'}
EXPECTED_VAL_IDS = ['VAL-'+str(n) for n in range(45, 55)]


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def classify(sample_id):
    if sample_id.startswith(PREFIX):
        return 'healthy_candidate', ''
    if sample_id in STANDARDS:
        return 'standard', STANDARDS[sample_id]
    return None


def identity(path):
    """Use the germline_seq sample folder when available; expose conflicts."""
    stem = next((path.name[:-len(s)] for s in SUFFIXES
                 if path.name.endswith(s)), None)
    if stem is None:
        return None
    folder = path.parent.parent.name if path.parent.name == 'germline_seq' else ''
    named, housed = classify(stem), classify(folder)
    if not named and not housed:
        return None
    sample_id = folder if housed else stem
    cohort, standard_id = housed or named
    # Matching an eligible folder does not authorize silently relabeling a file.
    conflict = bool(housed and stem != folder)
    return dict(sample_id=sample_id, cohort=cohort, standard_id=standard_id,
                identity_basis='folder_and_filename' if housed and not conflict
                else 'folder_name_conflict' if conflict else 'filename_only',
                filename_sample_id=stem, identity_conflict=conflict)


def inspect_file(path, root, identified):
    row = dict(identified, path=str(path), run_id=path.relative_to(root).parts[0],
               resolved_path='', size_bytes=None, modified_utc='',
               symlink=path.is_symlink(), physical_id='', readable=False,
               regular_nonempty=False, index_path='', index_readable=False,
               error='', header_checked=False)
    try:
        info = path.stat()
        row.update(resolved_path=str(path.resolve()), size_bytes=info.st_size,
                   modified_utc=datetime.fromtimestamp(info.st_mtime, timezone.utc).isoformat(),
                   physical_id=str(info.st_dev)+':'+str(info.st_ino),
                   readable=os.access(path, os.R_OK),
                   regular_nonempty=stat.S_ISREG(info.st_mode) and info.st_size > 0)
        for suffix in ('.tbi', '.csi'):
            index = Path(str(path)+suffix)
            if index.is_file():
                row['index_path'] = str(index)
                row['index_readable'] = os.access(index, os.R_OK)
                break
    except OSError as exc:
        row['error'] = str(exc)
    return row


def scan(root, progress=None):
    root = Path(root).expanduser().absolute()
    if not root.is_dir():
        raise ValueError('搜尋根目錄不存在或無法讀取：'+str(root))
    started = timestamp()
    files, errors, skipped_links = [], [], []
    directories = 0

    def onerror(exc):
        errors.append(dict(path=str(exc.filename or root), error=str(exc)))

    for current, dirs, names in os.walk(root, topdown=True, followlinks=False,
                                       onerror=onerror):
        directories += 1
        kept = []
        for name in sorted(dirs):
            path = Path(current)/name
            if path.is_symlink():
                skipped_links.append(str(path))
            else:
                kept.append(name)
        dirs[:] = kept
        for name in sorted(names):
            path = Path(current)/name
            identified = identity(path)
            if identified:
                row = inspect_file(path, root, identified)
                files.append(row)
                if row['error']:
                    errors.append(dict(path=str(path), error=row['error']))
        if progress and directories % 100 == 0:
            progress(directories, len(files))

    groups, physical = {}, {}
    for row in files:
        groups.setdefault(row['sample_id'], []).append(row)
        if row['regular_nonempty'] and row['physical_id']:
            physical.setdefault(row['physical_id'], set()).add(row['sample_id'])
    shared = {key: sorted(ids) for key, ids in physical.items() if len(ids) > 1}
    samples = []
    for sid, rows in sorted(groups.items()):
        usable = [r for r in rows if r['regular_nonempty'] and not r['identity_conflict']]
        physical_ids = {r['physical_id'] for r in usable}
        conflict = any(r['identity_conflict'] for r in rows)
        shared_ids = any(r['physical_id'] in shared for r in rows)
        selection = ('CHECK_IDENTITY' if conflict or shared_ids else
                     'NO_NONEMPTY_FILE' if not usable else
                     'MULTIPLE_FILES_NEED_REVIEW' if len(physical_ids) > 1 else
                     'SAME_FILE_MULTIPLE_PATHS' if len(usable) > 1 else 'ONE_FILE')
        samples.append(dict(sample_id=sid, cohort=rows[0]['cohort'],
                            standard_id=rows[0]['standard_id'], file_count=len(rows),
                            nonempty_matching_file_count=len(usable),
                            distinct_physical_files=len(physical_ids), selection=selection,
                            paths=[r['path'] for r in rows]))
        for row in rows:
            row['sample_selection'] = selection
            row['shared_file_across_ids'] = shared_ids
            row['ready_for_header_check'] = bool(row['regular_nonempty'] and row['readable']
                and row['index_readable'] and not row['identity_conflict'] and not shared_ids)
    health = [r for r in samples if r['cohort'] == 'healthy_candidate']
    controls = [r for r in samples if r['cohort'] == 'standard']
    counts = dict(healthy_sample_ids=len(health),
                  healthy_sample_ids_with_nonempty_matching_files=sum(
                      r['nonempty_matching_file_count'] > 0 for r in health),
                  healthy_gvcf_paths=sum(r['cohort']=='healthy_candidate' for r in files),
                  healthy_ids_with_multiple_files=sum(
                      r['distinct_physical_files'] > 1 for r in health),
                  standard_sample_ids=len(controls),
                  standard_biological_ids=len({r['standard_id'] for r in controls}),
                  standard_gvcf_paths=sum(r['cohort']=='standard' for r in files),
                  distinct_nonempty_physical_files=len(physical),
                  identity_conflicts=sum(r['identity_conflict'] for r in files),
                  physical_files_shared_across_ids=len(shared),
                  paths_ready_for_header_check=sum(r['ready_for_header_check'] for r in files),
                  errors=len(errors), skipped_directory_symlinks=len(skipped_links))
    return dict(schema='PRS_GVCF_INVENTORY_V1', started_at=started, finished_at=timestamp(),
                root=str(root), healthy_prefix=PREFIX,
                discovery_complete=not errors and not skipped_links,
                directories_scanned=directories, counts=counts,
                expected_val_ids_not_found=[s for s in EXPECTED_VAL_IDS
                    if s not in groups or not any(r['regular_nonempty'] and
                       not r['identity_conflict'] for r in groups[s])],
                samples=samples, files=files, errors=errors,
                skipped_directory_symlinks=skipped_links,
                physical_files_shared_across_ids=shared,
                scope=['Sample ID counts are not deduplicated biological people.',
                       '26T0 is the user-defined candidate rule, not confirmed disease-free status.',
                       'No gVCF content/header, genotype, genome coverage, ancestry or PRS checked.',
                       'No files selected automatically; no analysis or source changes.'])


def is_inside(path, root):
    try:
        Path(path).resolve().relative_to(Path(root).resolve())
        return True
    except ValueError:
        return False


def write_reports(report, output_dir):
    output_dir = Path(output_dir).expanduser().absolute()
    if is_inside(output_dir, report['root']):
        raise ValueError('報告請寫到既有 PRS work 目錄，不要寫進原始資料搜尋範圍。')
    output_dir.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix='inventory_'+datetime.now(timezone.utc).strftime(
        '%Y%m%dT%H%M%SZ')+'_', dir=output_dir))
    with (directory/'inventory.json').open('x', encoding='utf-8') as f:
        json.dump(report, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')
    for name, values, fields in (
        ('files.tsv', report['files'], ['cohort','sample_id','standard_id','run_id','path',
         'size_bytes','modified_utc','regular_nonempty','readable','index_path','index_readable',
         'identity_basis','filename_sample_id','identity_conflict','sample_selection',
         'symlink','resolved_path','physical_id','shared_file_across_ids','ready_for_header_check',
         'header_checked','error']),
        ('samples.tsv', report['samples'], ['cohort','sample_id','standard_id','file_count',
         'nonempty_matching_file_count','distinct_physical_files','selection','paths'])):
        with (directory/name).open('x', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for row in values:
                writer.writerow({**row, **({'paths':json.dumps(row['paths'],ensure_ascii=False)}
                                           if 'paths' in row else {})})
    return directory


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('/home/datalake_Raw/Novaseq'))
    parser.add_argument('--output-dir', type=Path, required=True,
                        help='在此既有 PRS work 目錄下新建報告子目錄，每次均不覆寫')
    args = parser.parse_args(argv)
    if is_inside(args.output_dir, args.root):
        parser.error('--output-dir 不能在原始資料搜尋範圍內')
    try:
        print('[PRS inventory] 唯讀盤點；不讀取 gVCF 內容、不執行分析。', flush=True)
        report = scan(args.root, lambda dirs, found: print(
            '[PRS inventory] directories='+str(dirs)+' matched_paths='+str(found), flush=True))
        directory = write_reports(report, args.output_dir)
    except (ValueError, OSError) as exc:
        print('[PRS inventory] ERROR: '+str(exc), file=sys.stderr)
        return 1
    c = report['counts']
    print('盤點狀態：'+('搜尋完成' if report['discovery_complete'] else
          '不完整；有讀取錯誤或未追蹤的目錄連結，數量僅為已找到的下限'))
    print('26T0 健康人候選檢體 ID：'+str(c['healthy_sample_ids']))
    print('其中有非空且名稱相符的檔案：'+str(c['healthy_sample_ids_with_nonempty_matching_files']))
    print('26T0 gVCF 路徑數：'+str(c['healthy_gvcf_paths']))
    print('同一健康人候選 ID 有不同實體檔案，待選來源：'+str(c['healthy_ids_with_multiple_files']))
    print('標準品檢體 ID／基準個體／gVCF 路徑：'+str(c['standard_sample_ids'])+'/'+
          str(c['standard_biological_ids'])+'/'+str(c['standard_gvcf_paths']))
    print('尚未找到的 VAL 標準品：'+(', '.join(report['expected_val_ids_not_found']) or '無'))
    print('名稱衝突／跨 ID 共用實體檔案：'+str(c['identity_conflicts'])+'/'+
          str(c['physical_files_shared_across_ids']))
    print('讀取錯誤／未追蹤目錄連結：'+str(c['errors'])+'/'+str(c['skipped_directory_symlinks']))
    print('報告：'+str(directory))
    print('以上是檢體 ID 去重數；個人去重、樣本身分、WGS 品質與 PRS 仍待核對。')
    return 0 if report['discovery_complete'] else 2


if __name__ == '__main__':
    raise SystemExit(main())

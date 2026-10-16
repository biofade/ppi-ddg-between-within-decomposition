# -*- coding: utf-8 -*-
"""scripts/v8_1/151_apply_author_metadata.py -- 用真实 metadata 填充首页/致谢/availability。

用法：
    1) 复制 FINAL_SUBMISSION/AUTHOR_METADATA.template.json 为
       FINAL_SUBMISSION/AUTHOR_METADATA.json 并逐字段填写真实信息
    2) 运行本脚本（可加 --check 只校验不写入）

纪律：
    * 任何字段为空 / 仍含 "TODO" / 含 "PLACEHOLDER" 的，**拒绝写入**并保持占位符，
      以免把无效信息印进终稿。
    * 本脚本绝不生成或猜测任何作者信息。
"""
import json
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".")
META = R / "FINAL_SUBMISSION/AUTHOR_METADATA.json"
TEX = R / "latex_v6/manuscript_v6.tex"
CHECK_ONLY = "--check" in sys.argv

if not META.exists():
    print("  未找到 AUTHOR_METADATA.json")
    print("  请先复制模板并填写：")
    print("    copy FINAL_SUBMISSION\\AUTHOR_METADATA.template.json "
          "FINAL_SUBMISSION\\AUTHOR_METADATA.json")
    print()
    print("  状态保持：UNRESOLVED — REQUIRES AUTHOR CONFIRMATION")
    sys.exit(2)

M = json.loads(META.read_text(encoding="utf-8"))


def bad(v):
    if v is None:
        return True
    s = str(v).strip()
    return (not s) or ("TODO" in s.upper()) or ("PLACEHOLDER" in s.upper())


def pending(v):
    return str(v).strip().upper().startswith("TODO")


problems = []
authors = M.get("authors") or []
aff_keys = M.get("affiliation_keys") or {}
if not authors:
    problems.append("authors 为空")
else:
    for i, a in enumerate(authors):
        for k in ("given", "family"):
            if bad(a.get(k)):
                problems.append(f"authors[{i}].{k} 未填")
        if not bad(a.get("orcid")) and not re.match(r"^\d{4}-\d{4}-\d{4}-[\dXx]{4}$",
                                                    str(a["orcid"]).strip()):
            problems.append(f"authors[{i}].orcid 格式不是 0000-0000-0000-0000")
        key = a.get("affiliation_key")
        if key is not None and bad(aff_keys.get(key)):
            problems.append(f"authors[{i}].affiliation_key='{key}' 未在 affiliation_keys 中定义")
    if len(authors) > 1 and not any(a.get("is_corresponding") for a in authors):
        problems.append("多作者但未指定通讯作者（is_corresponding=true）")

# 单一命名单位时不需要编号列表；否则必须有 affiliations
multi_aff = (not aff_keys) or len([k for k in aff_keys if not bad(aff_keys[k])]) > 1
affs = [v for v in (M.get("affiliations") or []) if not bad(v)]
if not aff_keys and not affs:
    problems.append("affiliations / affiliation_keys 均未填")

if bad(M.get("corresponding_email")):
    problems.append("corresponding_email 未填")
elif not re.match(r"^[^@\s]+@[^@\s]+\.[^@\s]+$", str(M["corresponding_email"]).strip()):
    problems.append("corresponding_email 格式不合法")

if bad(M.get("funding")):
    problems.append("funding 未填（无资助请写 'This work received no specific funding.'）")

repo = M.get("repository") or {}
route = str(repo.get("route", "")).strip()
url_pending = pending(repo.get("url", "")) or pending(repo.get("doi", "")) \
    or pending(repo.get("accession", ""))

# URL 仍是 TODO 时：按 'with_manuscript' 措辞写入（填了真实 URL 后重跑会自动升级措辞）
if url_pending and route == "public_repository":
    route = "with_manuscript"
    print("  注：repository.url 尚未提供（TODO），本次按 'with_manuscript' 措辞写入；")
    print("      你给出 GitHub URL 后重跑本脚本即可自动改为 public_repository 措辞。")

if route == "public_repository":
    if bad(repo.get("url")) and bad(repo.get("doi")) and bad(repo.get("accession")):
        problems.append("repository.route=public_repository 但 url/doi/accession 全未填")
elif route in ("with_manuscript", "on_request"):
    pass
else:
    problems.append("repository.route 非法（应为 public_repository / with_manuscript / "
                    "on_request）")

if problems:
    print("  拒绝写入 —— 以下字段尚未提供真实内容：")
    for p in problems:
        print("    -", p)
    print()
    print("  状态保持：UNRESOLVED — REQUIRES AUTHOR CONFIRMATION")
    sys.exit(1)

# ---------------- 生成 LaTeX 片段 ----------------
def tex_escape(s):
    return (str(s).replace("\\", r"\textbackslash{}")
            .replace("&", r"\&").replace("%", r"\%").replace("$", r"\$")
            .replace("#", r"\#").replace("_", r"\_")
            .replace("{", r"\{").replace("}", r"\}"))


# 姓名（不转义下划线以外字符，保持可读）
def nm(a):
    return f"{a['given']} {a['family']}"


# 单位：单一命名单位时不加编号；多单位时按序编号
if aff_keys and len([k for k in aff_keys if not bad(aff_keys[k])]) == 1:
    single_key = [k for k in aff_keys if not bad(aff_keys[k])][0]
    aff_of = {single_key: tex_escape(aff_keys[single_key])}
    multi = False
else:
    aff_of = {}
    multi = True
    for i, x in enumerate(affs):
        aff_of[f"__idx{i}"] = tex_escape(x)

author_parts = []
for a in authors:
    marks = []
    if multi:
        key = a.get("affiliation_key")
        if key and key in aff_keys:
            idx = [k for k in aff_keys if not bad(aff_keys[k])].index(key) + 1
            marks.append(str(idx))
        else:
            marks.append(str(a.get("affiliation_index", 1)))
    if len(authors) > 1 and a.get("is_corresponding"):
        marks.append("*")
    sup = ",".join(marks)
    author_parts.append(f"{nm(a)}$^{{{sup}}}$" if sup else nm(a))
authors_line = ", ".join(author_parts)

if multi:
    aff_lines = "\n".join(f"$^{{{i+1}}}$ {aff_of[f'__idx{i}']}" + r"\\"
                          for i in range(len(aff_of)))
else:
    aff_lines = list(aff_of.values())[0] + r"\\"

corr = [nm(a) for a in authors if a.get("is_corresponding")] or [nm(authors[0])]
corr_names = ", ".join(corr)

_orcids = [tex_escape(str(a["orcid"])) for a in authors if not bad(a.get("orcid"))]
thanks_body = aff_lines + "\n" + r"Corresponding author: " + tex_escape(corr_names) \
    + r" $\cdot$ Email: \texttt{" + tex_escape(M["corresponding_email"]) + "}"
if _orcids:
    thanks_body += "\n" + r"$\cdot$ ORCID: " + ", ".join(_orcids)

TITLE_BLOCK = ("\\author{" + authors_line + "\n"
               + r"\thanks{" + thanks_body + "\n}\n}\n")

FUNDING_TEX = tex_escape(M["funding"]) + " "

# ---------------- availability ----------------
if route == "public_repository":
    bits = []
    if not bad(repo.get("url")):
        bits.append(r"\url{" + str(repo["url"]).strip() + "}")
    if not bad(repo.get("doi")):
        bits.append("DOI: \\texttt{" + tex_escape(repo["doi"]) + "}")
    if not bad(repo.get("accession")):
        bits.append("accession: \\texttt{" + tex_escape(repo["accession"]) + "}")
    if not bad(repo.get("license")):
        bits.append("license: " + tex_escape(repo["license"]))
    AVAIL = ("All analyses use existing predictions and published external predictions; no new "
             "model was trained and no GPU was used. Frozen input hashes, the machine-readable "
             "source of every headline number and an automated consistency checker are "
             "available at " + "; ".join(bits) + ", and are also provided with the "
             "Supplementary Information.")
elif route == "with_manuscript":
    AVAIL = ("All analyses use existing predictions and published external predictions; no new "
             "model was trained and no GPU was used. Frozen input hashes, the machine-readable "
             "source of every headline number and an automated consistency checker are "
             "provided with the Supplementary Information and with the manuscript. A public "
             "repository will be linked here upon acceptance.")
else:
    AVAIL = ("All analyses use existing predictions and published external predictions; no new "
             "model was trained and no GPU was used. Frozen input hashes and the machine-readable "
             "source of every headline number are available from the corresponding author upon "
             "reasonable request.")

if CHECK_ONLY:
    print("  --check：校验通过，未写入")
    print("  作者:", authors_line)
    print("  通讯:", corr_names, "|", M["corresponding_email"])
    print("  route:", route)
    sys.exit(0)

s = TEX.read_text(encoding="utf-8")
orig = s

# 1) 作者块（lambda 避免 re.sub 解释替换串里的反斜杠）
#    占位块形如：
#      \author{Author Name(s)\thanks{Affiliation(s). Corresponding author: [NAME] $\cdot$
#      Email: [EMAIL] $\cdot$ ORCID: [ORCID\_PLACEHOLDER].}}
s, n_auth = re.subn(
    r"\\author\{Author Name\(s\)\\thanks\{.*?\}\}\n",
    lambda m: TITLE_BLOCK, s, count=1, flags=re.S)
if not n_auth:
    print("  警告：未匹配到作者占位块，作者信息未写入")
# 2) 致谢 funding
s = s.replace("[FUNDING\\_PLACEHOLDER] ", FUNDING_TEX)
# 3) data availability 段
s = re.sub(
    r"(\\section\*\{Data and code availability\}\n%={2,}\n).*?(\n\n%={2,}\n\\section\*\{Author contributions\})",
    lambda m: m.group(1) + AVAIL + m.group(2), s, count=1, flags=re.S)

if s == orig:
    print("  警告：未匹配到任何占位符，未写入任何内容")
    sys.exit(1)

TEX.write_text(s, encoding="utf-8")
print("  已写入作者/致谢/availability")
print("  作者:", authors_line)
print("  通讯:", corr_names, "|", M["corresponding_email"])
print("  funding:", M["funding"][:70])
print("  availability route:", route)
print()
print("  下一步：重建 PDF 并运行 validator")

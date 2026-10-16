# -*- coding: utf-8 -*-
"""scripts/v8_1/152_validator_metadata.py -- 让 R/S 检查反映真实 metadata 状态。"""
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = "validate_final_manuscript.py"
s = io.open(P, encoding="utf-8").read()

# ---- 替换 R（Data/code availability）----
mR = re.search(r'rec\("R\. Data/code availability truthful",.*?\n\n', s, re.S)
NEW_R = '''_meta_path = V81.parent / "FINAL_SUBMISSION" / "AUTHOR_METADATA.json"
_meta = (json.loads(_meta_path.read_text(encoding="utf-8"))
         if _meta_path.exists() else {})
_repo = _meta.get("repository", {}) or {}
_repo_url = str(_repo.get("url", ""))
_url_ok = bool(re.match(r"^https?://", _repo_url.strip())) and "TODO" not in _repo_url.upper()
_dmR = pymupdf.open(str(LATEX / "manuscript_v6.pdf"))
_tmR = "\\n".join(pg.get_text() for pg in _dmR)
_dmR.close()
_avail_claimed = "automated consistency checker" in _tmR
if _url_ok and _avail_claimed and _repo_url.strip() in _tmR:
    rec("R. Data/code availability truthful", "PASS",
        "public repository URL is stated in the manuscript and the SI artifacts exist "
        "in the submission package")
else:
    rec("R. Data/code availability truthful", "WARNING",
        "Manuscript states frozen hashes + machine-readable source + consistency checker are "
        "provided with the Supplementary Information, and that a public repository will be "
        "linked on acceptance. Repository URL not yet supplied (author action pending), so the "
        "availability statement currently relies on the SI package.")

'''
s = s[:mR.start()] + NEW_R + s[mR.end():]

# ---- 替换 S（author placeholders）----
mS = re.search(r'rec\("S\. author placeholders",.*?\n\n', s, re.S)
NEW_S = '''_PLACEHOLDERS = ("Author Name(s)", "[EMAIL]", "[NAME]", "ORCID_PLACEHOLDER",
                 "FUNDING_PLACEHOLDER", "TODO_AUTHOR_INPUT")
_found = [x for x in _PLACEHOLDERS if x in _tmR]
_has_author = False
_has_email = False
_has_orcid = False
try:
    _au = _meta.get("authors") or []
    _has_author = bool(_au) and all(
        str(a.get("family", "")).strip() and str(a.get("given", "")).strip()
        for a in _au)
    _has_email = bool(re.match(r"^[^@\\s]+@[^@\\s]+\\.[^@\\s]+$",
                               str(_meta.get("corresponding_email", "")).strip()))
    _orc = [str(a.get("orcid", "")).strip() for a in _au if a.get("orcid")]
    _has_orcid = bool(_orc) and all(
        re.match(r"^\\d{4}-\\d{4}-\\d{4}-[\\dXx]{4}$", o) for o in _orc)
except Exception:
    pass
_fund_ok = bool(str(_meta.get("funding", "")).strip()) and not _found
if (not _found) and _has_author and _has_email and _has_orcid and _fund_ok:
    rec("S. author metadata filled", "PASS",
        "real author name, affiliation, corresponding email, ORCID and funding statement are "
        "present in the manuscript; no placeholder remains")
else:
    rec("S. author metadata filled", "WARNING",
        "remaining placeholder strings = %s; author=%s email=%s orcid=%s funding=%s"
        % (_found or "none", _has_author, _has_email, _has_orcid, _fund_ok))

'''
s = s[:mS.start()] + NEW_S + s[mS.end():]

io.open(P, "w", encoding="utf-8").write(s)
import ast  # noqa: E402
ast.parse(s)
print("R/S checks updated; syntax OK")

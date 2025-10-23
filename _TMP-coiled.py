# ------------------------------------
# 
# scratch - fake run
# 
# debug coiled
# 
# ------------------------------------
import os
import sys
import json

import s3fs

try:
    from ew_workflows import scepter_helperFxns as shf
except Exception:
    shf = None


# --- 


def _coerce_value(val):
    """Try to coerce common literal strings to Python types.
    Falls back to the original string when no conversion applies.
    """
    if not isinstance(val, str):
        return val
    low = val.lower()
    if low in ("true", "yes", "on"):
        return True
    if low in ("false", "no", "off"):
        return False
    if low in ("none", "null"):
        return None
    # try integer/float
    try:
        if "." in val:
            return float(val)
        return int(val)
    except Exception:
        return val


def parse_dash_args(argv=None):
    """Parse a list of argv-style tokens into a dict.

    Rules:
    - --key value  -> {"key": value}
    - --key=value  -> {"key": value}
    - --flag        -> {"flag": True}
    - repeated keys -> collected into a list
    - positional args -> stored under "_positional" as a list
    - keys normalized: dashes -> underscores
    - values coerced by _coerce_value
    """
    if argv is None:
        argv = sys.argv[1:]
    parsed = {}
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok.startswith("--"):
            keyval = tok[2:]
            if "=" in keyval:
                key, raw = keyval.split("=", 1)
                val = _coerce_value(raw)
            else:
                key = keyval
                # lookahead for a value token
                if i + 1 < len(argv) and not argv[i + 1].startswith("--"):
                    val = _coerce_value(argv[i + 1])
                    i += 1
                else:
                    val = True
            key = key.replace("-", "_")
            # accumulate repeated keys into lists
            if key in parsed:
                if isinstance(parsed[key], list):
                    parsed[key].append(val)
                else:
                    parsed[key] = [parsed[key], val]
            else:
                parsed[key] = val
        else:
            parsed.setdefault("_positional", []).append(_coerce_value(tok))
        i += 1
    return parsed

print("Attempting to parse argv...")
# --- create an empty file TROUBLESHOOT ---
fs = s3fs.S3FileSystem()
with fs.open('s3://carbonplan-carbon-removal/ew-workflows-data/_flag.txt', 'w') as f:
    pass  # creates an empty file

parsed = parse_dash_args(sys.argv)
# print in JSON so it's easy to consume in tests or logs
print(json.dumps(parsed, indent=2, ensure_ascii=False))

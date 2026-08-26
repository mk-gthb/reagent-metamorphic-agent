#!/usr/bin/env python3
"""re:AGENT smoke test — confirms Proto + Modal + ESMFold work end to end,
and prints the API details needed to finalize T0/R9. Run: python smoke_test.py"""
import sys, inspect
import proto_tools

# 1) discover ESMFold entry points (robust to exact naming)
names = [x for x in dir(proto_tools) if "fold" in x.lower()]
print("== fold-related names in proto_tools ==\n", names)

run_fn = next((getattr(proto_tools, n) for n in names
               if n.lower().startswith("run_") and "fold" in n.lower()), None)
if run_fn is None:
    print("!! no run_*fold* function found — paste the list above and stop here.")
    sys.exit(0)
print("\n== signature of", run_fn.__name__, "==\n", inspect.signature(run_fn))

InputCls = next((getattr(proto_tools, n) for n in dir(proto_tools)
                 if "fold" in n.lower() and n.endswith("Input")), None)
ConfigCls = next((getattr(proto_tools, n) for n in dir(proto_tools)
                  if "fold" in n.lower() and n.endswith("Config")), None)
print("Input class:", InputCls, "| Config class:", ConfigCls)

# 2) run one prediction
seq = "MKTAYLLIGLLAIAAFSPQVLA"
try:
    out = run_fn(InputCls(sequences=[seq]), ConfigCls(device="modal"))
except Exception as e:
    print("\n!! run failed:", repr(e))
    print(">> paste this error and stop here.")
    sys.exit(0)

# 3) show where pLDDT + structure live (this is what we need to finalize the scripts)
print("\n== output structure ==")
print("type(out):", type(out))
print("dir(out):", [a for a in dir(out) if not a.startswith('_')])
res0 = out.results[0]
print("type(results[0]):", type(res0))
print("fields of results[0]:", [a for a in dir(res0) if not a.startswith('_')])
for attr in ("plddt", "mean_plddt", "confidence", "pdb", "pdb_string", "structure", "cif"):
    if hasattr(res0, attr):
        print(f"  results[0].{attr} -> {type(getattr(res0, attr))} "
              f"{str(getattr(res0, attr))[:80]}")
print("\nDONE — paste everything above.")

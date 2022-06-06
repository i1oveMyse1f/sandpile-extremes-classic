import os
import subprocess

to_compile_files = ["distribution", "conditional_probability_L_Eta_A_Kbins", "roc_curve_L_Eta_A", "event_rate"]

for f in to_compile_files:
    inp = os.path.join(".", "src", f + ".cpp")
    out = os.path.join(".", "bin", f)
    os.system(' '.join(["g++", "-std=c++17", "-o " + out, inp]))
    #subprocess.run(["g++", "-std=c++17", "-o " + out, inp])

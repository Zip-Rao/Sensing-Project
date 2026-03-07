import sys
import os

print("Current directory:", os.getcwd())
print("\nPython path before:")
for p in sys.path:
    print(" ", p)

# Simulate what the notebook does
current_dir = os.getcwd()
parent_dir = os.path.dirname(current_dir)
src_path = os.path.join(parent_dir, 'src')
print("\nTrying to append src_path:", src_path)
sys.path.append(src_path)

print("\nPython path after:")
for p in sys.path:
    print(" ", p)

print("\nTrying to import from src.qubit...")
try:
    from src.qubit import TransmonQubit
    print("SUCCESS: Imported TransmonQubit")
except ImportError as e:
    print("FAILED:", e)
    import traceback
    traceback.print_exc()

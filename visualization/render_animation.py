from paraview.simple import *
from tqdm import tqdm
import os
import sys

# === CONFIGURATION ===
STATE_FILE = "shield_scene.pvsm"
OUTPUT_DIR = "frames"
RESOLUTION = [2560, 1600]
FPS = 120
TOTAL_DURATION = 10.0  # Seconds of video
# =====================

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print("Loading state and re-routing pipeline...")
    LoadState(STATE_FILE)
    view = GetActiveView()

    # 1. Find the Reader and the Filters that depend on it
    reader = None
    sources = GetSources()

    for key, proxy in sources.items():
        name_str = key[0] if isinstance(key, tuple) else key
        if hasattr(proxy, 'FileName'):
            reader = proxy
            break

    if not reader:
        print("Error: No data reader found in state file.")
        sys.exit(1)

    # 2. Create the Interpolator
    interpolator = TemporalInterpolator(Input=reader)

    # 3. Re-route the pipeline: 
    # Find any filter whose Input is the reader, and change it to the interpolator
    for key, proxy in sources.items():
        if hasattr(proxy, 'Input'):
            if proxy.Input == reader:
                proxy.Input = interpolator
                proxy.UpdatePipeline()
                print(f"Re-routed: {key[0] if isinstance(key, tuple) else key} -> Interpolator")

    # 4. Final View Setup
    Hide(reader, view)
    # We don't necessarily call Show(interpolator) because the downstream 
    # filters (Tubes/Edges) are already shown and now connected to it.

    # 5. Timing Logic
    t_steps = reader.TimestepValues
    t_start = t_steps[0]
    t_end = t_steps[-1]
    num_frames = int(TOTAL_DURATION * FPS)

    print(f"Rendering {num_frames} frames at {RESOLUTION[0]}x{RESOLUTION[1]}...")

    # 6. Render Loop
    for i in tqdm(range(num_frames), desc="Rendering", unit="frame"):
        # Map loop index to simulation time
        current_time = t_start + (i / float(num_frames - 1)) * (t_end - t_start)
        
        view.ViewTime = current_time
        
        # Explicitly update pipeline to prevent 'No temporal data' error
        interpolator.UpdatePipeline(current_time)
        
        frame_path = os.path.join(OUTPUT_DIR, f"frame.{i:04d}.png")
        SaveScreenshot(frame_path, view, 
                    ImageResolution=RESOLUTION, 
                    CompressionLevel=1)

    print(f"TODO: convert frames to video")
          
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nRender halted by user.")
    except Exception as e:
        print(f"\nError during render: {e}")
    finally:
        print(f"\nProcess complete. Frames in {OUTPUT_DIR}/")
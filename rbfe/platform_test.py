import os

try:
    import openmm as mm          # OpenMM 8+
except ImportError:
    from simtk import openmm as mm   # OpenMM <8


def get_openmm_platform():
    """
    Return (platform, platform_properties) for Context().
    - Prefer GPU (CUDA or OpenCL) if available.
    - Otherwise use CPU with multithreading based on Slurm/OMP.
    """

    # Number of threads for CPU fallback
    n_threads = int(
        os.environ.get("SLURM_CPUS_PER_TASK", os.environ.get("OMP_NUM_THREADS", "1"))
    )

    # Check if Slurm actually gave us a GPU
    cuda_visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    has_slurm_gpu = cuda_visible not in ("", "none", "-1")

    # Try GPU platforms first if Slurm says we have a GPU
    if has_slurm_gpu:
        for name in ("CUDA", "OpenCL"):
            try:
                platform = mm.Platform.getPlatformByName(name)
                print(f"Using GPU platform: {name} (CUDA_VISIBLE_DEVICES={cuda_visible})")
                # No special properties needed for most GPU runs
                return platform, {}
            except Exception:
                pass  # try next GPU platform

    # Fallback: CPU with multithreading
    cpu = mm.Platform.getPlatformByName("CPU")
    cpu.setPropertyDefaultValue("Threads", str(n_threads))
    print(f"Using CPU platform with {n_threads} threads")
    return cpu, {"Threads": str(n_threads)}

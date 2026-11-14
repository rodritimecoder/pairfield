# Reproducing the ARCHITECTURE.md numbers

Standalone probes, no dependencies. Each prints the tables quoted in ARCHITECTURE.md.

    rustc -O research/verify_metrics.rs -o /tmp/vm && /tmp/vm
    rustc -O research/verify_itd_and_reflections.rs -o /tmp/vi && /tmp/vi

`verify_itd_and_reflections.rs` includes the image-source room model that produced the
reflection tap table (section 5).

#!/usr/bin/env python3
"""
Wrapper CLI that applies an external Boltz2 runtime patch before delegating
to the official Boltz CLI implementation.
"""

from __future__ import annotations

import sys

from boltz2_external_patch import apply_runtime_patch, load_config_from_env


def main() -> int:
    cfg = load_config_from_env()
    if cfg.enabled:
        apply_runtime_patch(cfg)

    from boltz.main import cli as boltz_cli

    # Delegate all args exactly as provided (e.g. "predict ...").
    boltz_cli.main(args=sys.argv[1:], prog_name="boltz", standalone_mode=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


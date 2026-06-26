from __future__ import annotations

import argparse
import runpy
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="oracle", description="ORACLE workflow CLI")
    sub = parser.add_subparsers(dest="command")

    init = sub.add_parser("init", help="Create an ORACLE project workspace")
    init.add_argument("workdir", type=Path)

    sub.add_parser("merlino", help="Delegate to the current Merlino CLI during migration")
    return parser


def main(argv: list[str] | None = None) -> int:
    args, remainder = build_parser().parse_known_args(argv)
    if args.command == "init":
        from oracle_core.workspace import ensure_workspace

        ensure_workspace(args.workdir)
        print(f"Created ORACLE workspace: {args.workdir}")
        return 0
    if args.command == "merlino":
        sys.argv = ["merlino", *remainder]
        runpy.run_module("merlino", run_name="__main__")
        return 0
    build_parser().print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


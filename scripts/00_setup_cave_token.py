"""Set up and verify the CAVE auth token (one time per machine). See docs/01_DATA_ACCESS.md §1."""
from __future__ import annotations

import sys


def main() -> int:
    try:
        from caveclient import CAVEclient
    except ImportError:
        print("caveclient not installed. Run: uv sync --extra dev")
        return 1

    print("MICrONS / CAVE auth setup")
    print("-" * 40)
    print("1. Accept the Terms of Service (once), with a Google-associated email:")
    print("   https://global.daf-apis.com/sticky_auth/api/v1/tos/2/accept")
    print("2. This will open a browser to mint a token.\n")

    client = CAVEclient()
    try:
        client.auth.setup_token(make_new=True)  # prints/open a URL to log in and get a token
    except Exception as e:  # noqa: BLE001
        print(f"setup_token note: {e}")

    token = input("\nPaste your token here (or leave blank if already saved): ").strip()
    if token:
        client.auth.save_token(token=token, overwrite=True)

    # Verify.
    from inhibitome.config import CFG

    c = CAVEclient(CFG.datastack)
    c.version = CFG.materialization_version
    versions = sorted(c.materialize.get_versions())
    print(f"\n✅ Authenticated to '{CFG.datastack}'. Pinned version {CFG.materialization_version}.")
    print(f"   Available materialization versions: {versions}")
    if CFG.materialization_version not in versions:
        print(f"⚠️  Pinned version {CFG.materialization_version} is not currently live "
              f"(may be expired). Pick one of {versions} in config/pilot.yaml.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

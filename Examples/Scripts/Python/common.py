from pathlib import Path


def getOpenDataDetectorDirectory():
    """
    Returns path to ODD files

    Located here so that the sources location can be obtained. The ODD files are not necessarily installed.
    """
    return (
        Path(__file__).parent.parent.parent.parent / "thirdparty" / "OpenDataDetector"
    )

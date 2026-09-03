from pathlib import Path

import pytest

TESTS = Path(__file__).resolve().parent
DATA = TESTS / 'data'
GOLDEN = TESTS / 'golden'


@pytest.fixture(scope='session')
def data_dir() -> Path:
    return DATA


@pytest.fixture(scope='session')
def golden_dir() -> Path:
    return GOLDEN

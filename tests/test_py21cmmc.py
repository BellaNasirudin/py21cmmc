
from click.testing import CliRunner

from py21cmmc import drive_21cmMC
from py21cmmc.cli import main


def test_main():
    runner = CliRunner()
    result = runner.invoke(main, [])

    assert result.output == '()\n'
    assert result.exit_code == 0


def test_drive_21cmMC():
    assert drive_21cmMC([b'a', b'bc', b'abc']) == b'abc'

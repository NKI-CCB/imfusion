"""Tests for the imfusion.util.shell module."""

import pytest

from imfusion.util import shell


class TestWhich(object):
    """Tests for the which function."""

    def test_positive_example(self):
        """Tests example with existing executable."""

        assert shell.which('ls') is not None

    def test_negative_example(self):
        """Tests example with non-existent executable."""

        assert shell.which('non-existent-binary-path') is None


class TestCheckDependencies(object):
    """Tests for the check_dependencies function."""

    def test_positive_example(self):
        """Tests example with existing executable."""

        shell.check_dependencies(['ls'])

    def test_negative_example(self):
        """Tests example with non-existent executable."""

        with pytest.raises(ValueError):
            shell.check_dependencies(['non-existent-binary-path'])


class TestParseArguments(object):
    """Tests for the parse_arguments function."""

    def test_example(self):
        """Tests a call with all basic argument types."""

        arg_str = '--single 5 --flag --multiple test1 test2 -s a'
        args = shell.parse_arguments(arg_str)

        assert len(args) == 4
        assert args['--flag'] == ()
        assert args['--single'] == ('5', )
        assert args['--multiple'] == ('test1', 'test2')
        assert args['-s'] == ('a', )


class TestFlattenArguments(object):
    """Tests for the flatten_arguments function."""

    def test_example(self):
        """Tests a call with all basic argument types."""
        args = {
            '--flag': True,
            '--single': 5,
            '--multiple': ('test1', 'test2'),
            '-s': ('a', )}  # yapf: disable

        flat = shell.flatten_arguments(args)
        assert flat == [
            '--flag', '--multiple', 'test1', 'test2', '--single', '5', '-s',
            'a'
        ]


class TestRunCommand(object):
    """Tests for the run_command function."""

    def test_example(self, mocker):
        """Tests a simple call to ls."""

        mock = mocker.patch.object(shell.subprocess, 'check_call')
        shell.run_command(args=['ls', '-l'], stdout=None, stderr=None)
        mock.assert_called_once_with(['ls', '-l'], stdout=None, stderr=None)

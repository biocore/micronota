from unittest import TestCase, main, mock
import gzip

from click.testing import CliRunner

from micronota.commands._rfam import cli


class Tests(TestCase):
    @mock.patch('micronota.commands._rfam.filter_models', return_value=9)
    def test(self, mock_add_metadata):
        runner = CliRunner()
        with runner.isolated_filesystem():
            with open('rfam.cm', 'w')  as f:
                f.write('>abc\nATGC')
            for op in ['bacteria', 'archaea', 'eukarya', 'default']:
                result = runner.invoke(cli, ['rfam.cm', '--operation', op, 'outfile'])
            self.assertEqual(result.exit_code, 0)


if __name__ == '__main__':
    main()

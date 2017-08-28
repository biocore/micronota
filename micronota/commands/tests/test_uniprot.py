from unittest import TestCase, main, mock
import gzip

from click.testing import CliRunner

from micronota.commands._uniprot import cli


class Tests(TestCase):
    @mock.patch('micronota.commands._uniprot.add_metadata', return_value=9)
    def test(self, mock_add_metadata):
        runner = CliRunner()
        with runner.isolated_filesystem():
            with gzip.open('uniprot.gz', 'w') as f1, open('test', 'w')  as f2:
                f1.write(b'>abc\nATGC')
                f2.write('>efg\nATGC')
            result = runner.invoke(cli, ['uniprot.gz', 'test', 'outfile'])
            self.assertEqual(result.exit_code, 0)


if __name__ == '__main__':
    main()

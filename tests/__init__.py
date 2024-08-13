import unittest
import os


class TestSConsPipelineOutput(unittest.TestCase):
    def test_pipeline_output(self):
        output_dir = 'test_output'

        output_files = [
            'blast',
            'lineages.csv',
            'seqs.fasta',
            'seq_info.csv',
            'taxonomy.csv',
            ]

        expected = {}
        for date_dir in os.listdir(output_dir):
            if len(date_dir) == 8 and date_dir.isnumeric():
                out = f'{output_dir}/{date_dir}'
                self.assertTrue(os.path.isfile(out + '/SUCCESS'))
                expected.update({
                    f'{out}/dedup/1200bp/types': output_files,
                    f'{out}/dedup/1200bp/named/': output_files,
                    f'{out}/dedup/1200bp/named/filtered': [
                        'blast',
                        'taxonomy.csv',
                        'outliers.csv',
                        'unsorted.fasta'],
                    f'{out}/dedup/1200bp/named/filtered/trusted': output_files,
                    f'{out}/dedup/1200bp/named/filtered/types': output_files,
                    })

        # Verify the folder structure and contents
        for folder_path, contents in expected.items():
            self.assertTrue(
                os.path.isdir(folder_path),
                f"Folder {folder_path} does not exist")
            actual_contents = os.listdir(folder_path)
            for expected_file in contents:
                self.assertIn(
                    expected_file,
                    actual_contents,
                    f"{expected_file} not found in {folder_path}")
                file = os.path.join(folder_path, expected_file)
                self.assertTrue(os.path.isfile(file))
                self.assertTrue(os.stat(file).st_size > 0)  # not empty


if __name__ == '__main__':
    unittest.main()

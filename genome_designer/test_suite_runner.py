"""
Extends the nose runner to create a temporary filesystem during tests that is
subsequently torn down.
"""

import os
import shutil

from django_nose import NoseTestSuiteRunner

import settings


class TempFilesystemTestSuiteRunner(NoseTestSuiteRunner):
    """Subclasses the nose test runner in order to use a temp filesystem.
    """

    def setup_databases(self):
        self.setup_temp_filesystem()
        return super(TempFilesystemTestSuiteRunner, self).setup_databases()

    def teardown_databases(self, *args, **kwargs):
        self.teardown_temp_filesystem()
        return super(TempFilesystemTestSuiteRunner, self).teardown_databases(
                *args, **kwargs)

    def setup_temp_filesystem(self):
        settings.MEDIA_ROOT = settings.TEST_FILESYSTEM_DIR
        if os.path.exists(settings.MEDIA_ROOT):
            shutil.rmtree(settings.MEDIA_ROOT)
        os.mkdir(settings.MEDIA_ROOT)

    def teardown_temp_filesystem(self):
        assert settings.MEDIA_ROOT == settings.TEST_FILESYSTEM_DIR
        if os.path.exists(settings.MEDIA_ROOT):
            shutil.rmtree(settings.MEDIA_ROOT)

class CustomTestSuiteRunner(TempFilesystemTestSuiteRunner):
    """Our custom TestSuiteRunner.
    """

    def setup_test_environment(self, **kwargs):
        assert_no_orphaned_pyc_files('.')
        return super(TempFilesystemTestSuiteRunner, self).setup_test_environment()


def assert_no_orphaned_pyc_files(start_dir):
    """Walk from the start directory through python pkg dirs, looking for
    .pyc files with no corresponding .py file.

    This is to avoid the bug where we delete a module but the .pyc file stays
    and allows tests to pass on the local machine.

    Raises AssertionError if an orphaned .pyc file is found.
    """
    orphaned_files = []
    for (dirpath, dirnames, filenames) in os.walk(start_dir, topdown=True):
        # NOTE: Python docs say it is okay to modify dirnames in-place in order
        # to prune the directories that we are walking.
        if not '__init__.py' in filenames:
            # Don't recurse any further.
            while len(dirnames):
                dirnames.pop()
            # No need to check this dir.
            continue


        # Build a set of names ending in .py.
        py_set = set()
        for filename in filenames:
            if os.path.splitext(filename)[1] == '.py':
                py_set.add(os.path.splitext(filename)[0])

        # Check that all .pyc files have a corresponding .py file.
        for filename in filenames:
            if os.path.splitext(filename)[1] == '.pyc':
                if not os.path.splitext(filename)[0] in py_set:
                    full_orphan_path = os.path.join(dirpath, filename)
                    orphaned_files.append(full_orphan_path)

    if len(orphaned_files):
        raise AssertionError, (
                "The following files are orphaned .pyc files:\n\n%s\n\n" %
                        '\n'.join(orphaned_files) +
                "Perhaps you moved the file and meant to delete them?")

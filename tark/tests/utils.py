from django.conf import settings
from django.test.runner import DiscoverRunner

class ManagedModelTestRunner(DiscoverRunner):
    """
    Test runner that automatically makes all unmanaged models in your Django
    project managed for the duration of the test run, so that one doesn't need
    to execute the SQL manually to create them.
    """
    def setup_test_environment(self, *args, **kwargs):
        settings.TESTING_MODE = True
        settings.SEQ_HOST_URI = None
        seen_tables = []
        from django.apps import apps
        self.unmanaged_models = [m for m in apps.get_models()
                                 if not m._meta.managed]
        for m in self.unmanaged_models:
            print "Flipping status of table {}".format(m.__name__)

            if m.__name__ == 'Genenames':
                print "Setting unique=True for Genenames.external_id"
                field = m._meta.get_field('external_id')
                setattr(field, '_unique', True)

            if m._meta.db_table in seen_tables:
                print "Seen {} already, skipping".format(m._meta.db_table)
                continue
            seen_tables.append(m._meta.db_table)
            m._meta.managed = True
        super(ManagedModelTestRunner, self).setup_test_environment(*args,
                                                                   **kwargs)

    def teardown_test_environment(self, *args, **kwargs):
        super(ManagedModelTestRunner, self).teardown_test_environment(*args,
                                                                      **kwargs)
        # reset unmanaged models
        for m in self.unmanaged_models:
            m._meta.managed = False

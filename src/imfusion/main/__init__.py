class Command(object):
    """Base command class."""

    name = 'base'

    def configure(self, parser):
        """Configures the argument parser for the command."""
        raise NotImplementedError()

    def run(self, args):
        """Runs the command."""
        raise NotImplementedError()

    @classmethod
    def available_commands(cls):
        """Returns dict of available commands."""

        return {
            class_.name: class_()
            for class_ in cls._get_subclasses() if class_.name != 'base'
        }

    @classmethod
    def _get_subclasses(cls):
        for subclass in cls.__subclasses__():
            for sub_subclass in subclass._get_subclasses():
                yield sub_subclass
            yield subclass

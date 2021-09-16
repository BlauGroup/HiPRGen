from monty.json import MSONable

class DispatcherPayload(MSONable):
    """
    class for storing all the arguments required by the reaction
    filter dispatcher. We do this instead of passing arguments
    directly because it makes it easier to pass arguments through the
    MPI barrier.
    """

    def __init__(
            self,
            bucket_db_file,
            reaction_network_db_file,
            report_file,
            commit_frequency = 1000,
            checkpoint_interval = 10):

        self.bucket_db_file = bucket_db_file
        self.reaction_network_db_file = reaction_network_db_file
        self.report_file = report_file
        self.commit_frequency = commit_frequency
        self.checkpoint_interval = checkpoint_interval


class WorkerPayload(MSONable):
    """
    class for storing all the arguments required by the reaction
    filter dispatcher.
    """
    def __init__(
            self,
            bucket_db_file,
            reaction_decision_tree,
            params,
            logging_decision_tree):

        self.bucket_db_file = bucket_db_file
        self.reaction_decision_tree = reaction_decision_tree
        self.params = params
        self.logging_decision_tree = logging_decision_tree

class Flat_file_loop:

    '''Simple  base  class for  looping  over  plain  files and  fasta
    files.'''

    def __init__(self, filename=None, report_period=None, fileobj=None):

        if fileobj is None:
            self.file = open(filename)
        else:
            self.file = fileobj
        self.sequence_count=-1

        # How often to print how many sequences have been read in.
        self.report_period = report_period

    def __getitem__(self, index):
        if not self.read():
            raise IndexError, 'End of plain file'
        if self.report_period and \
           (not (self.sequence_count % self.report_period)):
            print self.sequence_count
        return self.name, self.sequence

    def tell(self):

        return self.file.tell()

class Fasta_file(Flat_file_loop):

    def __init__(self, filename=None,report_period=None,fileobj=None,
                 break_size=None):
        Flat_file_loop.__init__(self, filename, report_period,fileobj)
        self.next_name = self.file.readline().strip()
        self.break_size = break_size
        self.finished_p = False

    def readall(self):
        sequence = []
        total_length = 0
        for line in self.file:
            if line.startswith('>'):
                self.name = self.next_name
                self.next_name = line.strip()
                break
            sequence.append(line.strip())
            total_length += len(sequence[-1])
            if (self.break_size is not None) and \
               total_length > self.break_size:
                if not hasattr(self, 'name'):
                    self.name = self.next_name
                break
        else:
            if not sequence: # Reached end of file.
                self.next_name = None
                self.finished_p = True
                return None, None
            self.name = self.next_name
        self.sequence = ''.join(sequence)
        self.sequence_count += 1
        return self.name, self.sequence

    def read(self): return self.readall()[1]

    def tell(self):

        return self.file.tell() - len(self.file.line)

    def next(self):
        rv = self.readall()
        if rv == (None, None):
            raise StopIteration
        return rv

    def __iter__(self):
        return self

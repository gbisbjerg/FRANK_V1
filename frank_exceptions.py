class ReallyBadInput(Exception):
    '''Exception raised when the sample input file does not contain $DATA'''
    pass
  
class SigsegvError(Exception):
    '''Exception raised when GAMESS has an error reading/writing memory,
       meaning that the calculation needs to be rerun without changes'''
    pass

class SigtermError(Exception):
    '''Exception raised when the calculation was terminated, meaning that
       the job needs to be rerun without changes'''
    pass

class GamessMemoryError(Exception):
    '''Exception raised when the calculation was not allocated enough memory,
       meaning memory needs to be adjusted in the input card'''
    pass

class GamessError(Exception):
    '''Exception raised for unexpected GAMESS errors, meaning the input card 
       likely needs to be corrected'''
    pass

class SCFConvergenceError(Exception):
    '''Exception raised when a SCF calculation fails to converge, meaning that
       the job didn't truly work and needs to be rerun'''
    pass

class BugError(Exception):
    '''Exception raised when GAMESS doesn't know what the error was, meaning 
       that the job needs to be rerun'''

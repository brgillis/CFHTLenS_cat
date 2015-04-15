import numpy as np
import astropy.io.ascii as ascii

class columns_reader(object):

    def __init__(self):
        self._keys_ = []
        self._col_names_ = []
        self._offsets_ = []
        self._factors_ = []
        self._indices_from_keys_ = {}
        self._keys_from_indices_ = {}
        
        self._filename_ = None
        self._cols_changed_ = True
        self._read_in_data_ = []
        
        return
        
    def add(self, key, col_name, offset=0., factor=1.):
        
        i = len(self._col_names_)
        
        self._indices_from_keys_[key] = i
        self._keys_from_indices_[i] = key
        
        self._keys_.append(key)
        self._col_names_.append(col_name)
        self._offsets_.append(offset)
        self._factors_.append(factor)
        
        self._cols_changed_ = True
        
        return
        
    def read(self, filename):
        
        if((self._cols_changed_) or (self._filename_ != filename)):

            # Load in the table
            try:
                data_table = ascii.read(filename)
            except:
                print("ERROR: Table " + filename + " cannot be read.")
                return
            
            self._read_in_data_ = []
            for col_i in xrange(self.num_cols()):
                self._read_in_data_.append((np.array(data_table[self._col_names_[col_i]])+self._offsets_[col_i])*self._factors_[col_i])
            
            self._cols_changed_ = False
            self._filename_ = filename
            
        return self._read_in_data_
    
    def num_cols(self):
        return len(self._col_names_)
    
    def keys(self):
        return self._keys_
    
    def col_names(self):
        return self._col_names_
    
    def data(self):
        return self._read_in_data_
    
    def indices(self):
        return self._indices_from_keys_.values()
    
    def indices_from_keys(self):
        return self._indices_from_keys_()
    
    def keys_from_indices(self):
        return self._keys_from_indices_()
    
    def index(self, key):
        return self._indices_from_keys_[key]
    
    def key(self, index):
        return self._keys_from_indices_[index]
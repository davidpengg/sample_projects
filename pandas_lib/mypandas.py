import operator
from pdb import set_trace

VALID_TYPES = {str, bool, int, float}

class Series(object):
    '''
    Series:
    - holds values of a certain type, as well as the None value
    - types: strings, booleans, floating point, interger
    - raise exception for operations between Series of different types (unless explicitely mentioned otherwise)
    - raise exception if lengths of Series don't match
    - construction: given a list of itmes to contain (items should be type checked)
    - use the square bracket operator: integer for value in Series at that position, ???? boolean Series for a new Series with all rows values where boolean Series is true (the "index" of a value should be a dictionary mapping bc of this function request)

    https://docs.python.org/3/reference/datamodel.html
    https://docs.python.org/3/library/operator.html

    All operations should return a Series.

    Behaviors of Series types, from interpreting specs and a little bit of experimenting with pandas Series:

    op      ..  str bool    int float            
    ==      ->  y   y       y   y  
    &       ->  -   y       -?  -?
    |       ->  -   y       -?  -?
    ^       ->  -   y       -?  -?
    ~       ->  -   y       -?  -?
    
    +       ->  y   -?      y   y
    -       ->  -   -       y   y
    *       ->  -   -       y   y
    /       ->  -   -       y   y
    //      ->  -   -       y   y

    >       ->  -?  -       y   y
    <       ->  -?  -       y   y
    >=      ->  -?  -       y   y
    <=      ->  -?  -       y   y
    
    !=      ->  y   y       y   y
    
    '''
    
    def __init__(self, elements = None): 
        '''
        Create a Series object

        elements : array or dictionary of indices to elements. 
        Contains data stored in Series. Included option to use 
        dictionary in case of non continuous indices from
        future Series slicing. Ex. ['a', 'b', 'c', 'd'] or 
        {0: 'a', 1 : 'b', 2 : 'c', 3 : 'd'}
        '''
        
        if not elements:
            raise Exception("Series cannot be empty (for now).")
        
        self.len = len(elements)

        input_type = type(elements)
        if input_type not in [dict, list]:
            raise TypeError("Series input must be list or dict.")

        NoneType = type(None)
        self.type = NoneType

        # type checking
        # runs in linear time
        for e in elements:
            if input_type == list:
                etype = type(e)
            else: #input_type == dict
                etype = type(elements[e])
        
            if self.type == NoneType and etype != NoneType:
                self.type = etype
            elif etype != NoneType and etype != self.type:
                raise TypeError("Series only holds values of a certain type (and the None value)")

        if self.type == NoneType:
            raise TypeError("Series must contain values that are not None")
        if self.type not in VALID_TYPES:
            raise TypeError("Series must contain only one of str, bool, int, float (including None value)")
        
        # create the actual data struct
        # runs in linear time
        if input_type == list:
            self.data = {i : elements[i] for i in range(self.len)}    
        else: # index == dict:
            # set_trace()
            index = list(elements.keys())
            self.data = {k : elements[k] for k in elements}

    def __getitem__(self, key):
        '''
        Gets element of Series

        key : int or Series. If int, gets element at 
        the index specified. If Series of booleans,
        gets elements at indices where boolean is True.
        '''
        if type(key) == int:
            return self.data[key]

        elif type(key) == Series:
            if key.type != bool:
                raise KeyError("Accessing Series must be done with index or with Series of type bool.")
            elif key.len != self.len:
                raise IndexError("Accessing Series with a Series of type bool must be of same length.")
            
            new_self = {}

            # indexing using bool Series
            # runs in linear time
            for e in self.data:
                if key[e]:
                    new_self[e] = self.data[e]
            
            return Series(new_self)
        else:
            raise KeyError("Not a valid key.")

    def __same_type(self, other):
        if self.type != other.type:
            raise ValueError("Operations between Series must be same types.")

    def __same_len(self, other):
        if self.len != other.len:
            raise ValueError("Operations between Series must be same length.")
    
    # see Series __init__ description on speeding up syntactic sugar Series
    def __sugar_other(self, other):
        '''
        Turns single element into Series of necessary 
        length. Syntactic sugar

        other : single value to be turned into Series
        of same length as self. If not a single value,
        returns the same unmodified Series.
        '''
        if type(self) != Series:
            raise TypeError("Left value must be a Series object.")
        
        # syntactic sugaring
        if type(other) != Series:
            if type(other) not in VALID_TYPES:
                raise TypeError("Right value must be one of str, bool, int, float, or Series object.")
            sugar_other = [other] * self.len
            other = Series(sugar_other)
        
        return other

    def __validation(self, other):
        '''
        Checks whether types and lengths of Series
        to be compared are the same.
        '''
        other = self.__sugar_other(other)
        self.__same_type(other)
        self.__same_len(other)

        return other

    # overwrite class operation functions

    def __repr__(self):
        '''
        Format Series for printing.
        '''

        s = ''
        for index in sorted(self.data):
            s += f"{index}\t{self.data[index]}\n"
        
        # if self.name != '':
        #     s += f"name: {self.name}, "
        
        s += f"type: {self.type.__name__}"
        return s

    def __class_op(self, other, op):
        '''
        Generic Series operation.
        '''

        res = {e: op(self.data[e], other.data[e]) for e in self.data}
        # set_trace()
        # for i, e in enumerate(self.data):
        #     res[i] = op(self.data[e], other.data[e]) \
        #         if (self.data[i] is not None and \
        #             other.data[i] is not None) \
        #         else None

        return res

    def __eq__(self, other):
        '''
        Series == operation.
        '''
        other = self.__validation(other)

        equal = self.__class_op(other, operator.eq)

        return Series(equal)

    def __ne__(self, other):
        '''
        Series != operation.
        '''
        other = self.__validation(other)

        noteq = self.__class_op(other, operator.ne)

        return Series(noteq)

    def __and__(self, other):
        '''
        Series & operation.
        '''
        other = self.__validation(other)

        if self.type != bool:
            raise TypeError("The '&' operator is not supported for this type.")
        else:
            anding = self.__class_op(other, operator.and_)
        
        return Series(anding)

    def __or__(self, other):
        '''
        Series | operation.
        '''
        other = self.__validation(other)

        if self.type != bool:
            raise TypeError("The '|' operator is not supported for this type.")
        else:
            oring = self.__class_op(other, operator.or_)

        return Series(oring)
    
    def __xor__(self, other):
        '''
        Series ^ operation.
        '''
        other = self.__validation(other)

        if self.type != bool:
            raise TypeError("The '^' operator is not supported for this type.")
        else:
            xoring = self.__class_op(other, operator.xor)

        return Series(xoring)

    def __invert__(self):
        '''
        Series ~ operation.
        '''

        if self.type != bool:
            raise TypeError("The '~' operator is not supported for this type.")
        else:
            inverted = [True] * self.len
            for i, e in enumerate(self.data):
                inverted[i] = not self.data[e]
        
        return Series(inverted)

    def __add__(self, other):
        '''
        Series + operation.
        '''

        other = self.__validation(other)

        if self.type == bool:
            added = self.__class_op(other, operator.or_)
        else:
            added = self.__class_op(other, operator.add)
        
        return Series(added)

    # beginning of int, float type only
    def __sub__(self, other):
        '''
        Series - operation.
        '''

        other = self.__validation(other)

        if self.type not in [int, float]:
            raise TypeError("The '-' operator is not supported for this type.")
        else:
            subtracted = self.__class_op(other, operator.ub)
        
        return Series(subtracted)
    
    def __mul__(self, other):
        '''
        Series * operation.
        '''
        other = self.__validation(other)

        if self.type not in [int, float]:
            raise TypeError("The '*' operator is not supported for this type.")
        else:
            mult = self.__class_op(other, operator.mul)
        
        return Series(mult)
    
    # def __div__(self, other):

        # div = self.__class_op(other, operator.div)

    #     return Series(div)

    def __floordiv__(self, other):
        '''
        Series // operation.
        '''
        other = self.__validation(other)

        if self.type not in [int, float]:
            raise TypeError("The '//' operator is not supported for this type.")
        else:
            div = self.__class_op(other, operator.floordiv)
        
        return Series(div)
        
    def __truediv__(self, other):
        other = self.__validation(other)

        if self.type not in [int, float]:
            raise TypeError("The '/' operator is not supported for this type.")
        else:
            div = self.__class_op(other, operator.truediv)
        
        return Series(div)

    def __lt__(self, other):
        '''
        Series < operation.
        '''

        other = self.__validation(other)

        if self.type not in [int, float]:
            raise TypeError("The '<' operator is not supported for this type.")
        else:
            comp = self.__class_op(other, operator.lt)
        
        return Series(comp)

    def __le__(self, other):
        '''
        Series <= operation.
        '''
        
        other = self.__validation(other)

        if self.type not in [int, float]:
            raise TypeError("The '<=' operator is not supported for this type.")
        else:
            comp = self.__class_op(other, operator.le)
        
        return Series(comp)
    
    def __gt__(self, other):
        '''
        Series > operation.
        '''

        other = self.__validation(other)

        if self.type not in [int, float]:
            raise TypeError("The '>' operator is not supported for this type.")
        else:
            comp = self.__class_op(other, operator.gt)
        
        return Series(comp)
    
    def __ge__(self, other):
        '''
        Series >= operation.
        '''

        other = self.__validation(other)

        if self.type not in [int, float]:
            raise TypeError("The '>=' operator is not supported for this type.")
        else:
            comp = self.__class_op(other, operator.ge)
        
        return Series(comp)

class DataFrame(object):
    '''
    DataFrame
    
    Two-dimensional, size-mutable, potentially heterogeneous tabular data.
    A dictionary of Series.

    All operations should return a DataFrame.

    '''
    
    def __init__(self, elements = None):
        '''
        Creates DataFrame

        elements : dictionary of Series, where key is
        the name of the column, value is the Series object
        '''

        if not elements:
            raise Exception("DataFrame cannot be empty (for now).")

        n = -1
        if type(elements) == dict:
            for s in elements:
                if type(elements[s]) != Series:
                    raise TypeError("DataFrame should contain Series")
                
                if n == -1:
                    # first Series: all other Series must have same length, same indices
                    n = elements[s].len
                    self.index = elements[s].data.keys()
                else:
                    if n != elements[s].len:
                        raise Exception("Series in a DataFrame must be of same length.")  
                    
                    if elements[s].data.keys() != self.index:
                        raise Exception("Series in a DataFrame must have the same indices.")
                # elements[s].name = s

            self.columns = list(elements.keys())
            
            self.table = elements 

            self.index = list(self.index)
        else:
            raise ValueError("DataFrame constructor takes only dictionary of Series")

    def __getitem__(self, key):
        '''
        Gets subset of DataFrame

        key : str or Series of booleans. If str, 
        return Series of that column name. If Series,
        return DataFrame of rows at indices where the
        booelan is True.
        '''
        if type(key) == str:
            return self.table[key]

        elif type(key) == Series and key.type == bool and key.len == len(self.index):

            # format: {Series1_name: {0: e0, 1: e1, ...}, Series2_name: {0 : f0, 1: f1, ...}, ...}
            new_elements_dict = {}
            for s in self.columns:
                new_elements_dict[s] = {}

            for i in self.index:
                if key[i]:
                    for s in new_elements_dict:
                        new_elements_dict[s][i] = self.table[s][i]
            
            # format: {Series1_name: Series({0: e0, 1: e1, ...}), Series2_name: Series({0: f0, 1: f1, ...}), ...}
            new_elements_Series = {s: Series(new_elements_dict[s]) for s in new_elements_dict}

            return DataFrame(new_elements_Series)

        else:
            raise KeyError("Not a valid key.")

    def __repr__(self):
        '''
        Format DataFrame for printing.
        '''

        s = '\t'
        s += '\t'.join(self.columns) + '\n'
        
        for idx in self.index:
            s += f"{idx}\t"

            for key in self.columns[:-1]:
                s += f"{self.table[key][idx]}\t"
            s += f"{self.table[self.columns[-1]][idx]}\n"
        s = s[:-1]

        return s
    
def read_csv(path, separator=','):
    '''
    Reads .csv into a new DataFrame
    
    path : str, path
    '''

    # names = {}

    with open(path) as f:
        lines = f.readlines()
        lines = [line.split(separator) for line in lines]
    
    n = len(lines[0])
    names = lines[0]
    data = {}

    for line in lines[1:]:
        for i in range(n):
            if names[i] not in data:
                data[names[i]] = [line[i]]
            else:
                data[names[i]].append(line[i])
    
    fin_data = {name : Series(data[name]) for name in data}
    return DataFrame(fin_data)

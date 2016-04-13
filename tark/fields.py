from django.db import models
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

ALPHABETS = {
             'gene': IUPAC.ambiguous_dna,
             'transcript': IUPAC.ambiguous_dna,
             'exon': IUPAC.ambiguous_dna,
             'translation': IUPAC.extended_protein
    }

class ChecksumField(models.CharField):
    description = "Allow retrieval of binary mysql fields"
    __metaclass__ = models.SubfieldBase

    def __repr__(self):
        return ''.join( [ "%02X" % ord( x ) for x in super(ChecksumField, self).__repr__ ] ).strip()
    
    def db_type(self, connection):
        return 'binary(20)'
    
    def from_db_value(self, value, expression, connection, context):
        if value is None:
            return value
        return ''.join( [ "%02X" % ord( x ) for x in value ] ).strip()
    
    def value_to_string(self, obj):
        value = self._get_val_from_obj(obj)
        return self.get_prep_value(value)

    def get_prep_value(self, value):
        if not isinstance(value, str):
            value = value.seq_checksum
        bytes = []
        for i in range(0, len(value), 2):
            bytes.append( chr( int (value[i:i+2], 16 ) ) )

        return ''.join( bytes )

class SequenceField(models.TextField):
    @classmethod
    def alphabet_type(cls, feature_type):
        return ALPHABETS.get(feature_type, IUPAC.Alphabet)
    
    def from_db_value(self, value, expression, connection, context):
        if value is None:
            return Seq()
        
        return Seq(value)
    
    def value_to_string(self, obj):
        value = self._get_val_from_obj(obj)
        return self.get_prep_value(value)

    def to_python(self, value):
        return str(value)
    
    def get_prep_value(self, value):
        if value is None:
            return None

        return str(value)
    
    def get_prep_lookup(self, lookup_type, value):
        if value is None:
            return None
        
        if lookup_type == 'isnull' and value in (True, False):
            return value
        
        if isinstance(value, (list, set)):
            return [ self.item_field_type.to_python(x) for x in value ]
        
        return self.item_field_type.to_python(value)

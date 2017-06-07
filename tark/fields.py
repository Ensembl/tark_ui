from django.db import models
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pprint
from _ast import alias
from django.db.models.fields.related_descriptors import ForwardManyToOneDescriptor

ALPHABETS = {
             'gene': IUPAC.ambiguous_dna,
             'transcript': IUPAC.ambiguous_dna,
             'exon': IUPAC.ambiguous_dna,
             'translation': IUPAC.extended_protein
    }

class ChecksumField(models.CharField):
    description = "Allow retrieval of binary mysql fields"
#    __metaclass__ = models.SubfieldBase

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
        if not isinstance(value, (str, unicode, buffer)):
            value = value.seq_checksum
        bytes = []

        for i in range(0, len(value), 2):
            bytes.append( chr( int (value[i:i+2], 16 ) ) )

        return ''.join( bytes )

class GeneSetField(models.CharField):
    description = "Allows retrieval of abstract field type containing a grouping of genes"
    
    def from_db_value(self, value, expression, connection, context):
        if value is None:
            return value
        decoded_genes = []
        genes = value.split('#!#!#')

        for gene in genes:
            pieces = gene.split(',', 3)
            if len(pieces) >= 3:
                pieces[3] = ''.join( [ "%02X" % ord( x ) for x in pieces[3] ] ).strip()
                
                decoded_genes.append(pieces)
        
#        pprint.pprint(decoded_genes)
        
        return decoded_genes
 
    def to_python(self, value):
        return value
    
class SequenceField(models.TextField):
    @classmethod
    def alphabet_type(cls, feature_type):
        return ALPHABETS.get(feature_type, IUPAC.Alphabet)
    
    def from_db_value(self, value, expression, connection, context):
        if value is None:
            return Seq('')
        
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

class HGNCForwardManyToOneDescription(ForwardManyToOneDescriptor):
    """
    We need to address what occurs when a lookup in the gene_name table
    fails because the HGNC column is NULL. There might be a more eligant
    way to handle this in django, I just haven't found it. So we're
    overloading the relationship manager such that when a __get__ is
    done, if the related column in gene_name fails we return none
    rather than letting the exception bubble up.
    """
    def __get__(self, instance, cls=None):
        
        try:
            rel_obj = super(HGNCForwardManyToOneDescription, self).__get__(instance, cls)
        except Exception as e:
            return None

        return rel_obj
    
class HGNCField(models.ForeignKey):
    """
    Derived class from ForeignKey, adding the restriction when following the
    froreign key to HGNC names, we only want the primary name from source type
    HGNC. Overloading get_extra_descriptor_filter just returns extra join
    conditions when looking up in the gene_names table.
    """
    requires_unique_target = False

    def get_extra_descriptor_filter(self, instance):
        return {'primary_id': 1, 'source': 'HGNC'}

    def contribute_to_class(self, cls, name, private_only=False, **kwargs):
        '''
        Override the relationship manager used for our HGNCField, see HGNCForwardManyToOneDescription
        above for more details.
        '''
        super(HGNCField, self).contribute_to_class(cls, name, private_only, **kwargs)

        setattr(cls, self.name, HGNCForwardManyToOneDescription(self))

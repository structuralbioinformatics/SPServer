from xml.dom.minidom import parse
'''
ATENTION: In xml.dom.minidom library you have two types of objects, Node and NodeList. NodeList are iterable
          lists. Nodes can be either Element and Text ones (among others like Attribute that are not used here). 
          Text nodes contains the text of the edges of the tree (and can be also in the core, as we will see in 
          the example).Element nodes are the nodes defined by a tag.
          Example:
          <author>
             <name>M.A.M.L</name>
          </author>
          The M.A.M.L is a text node which shoud be accesed using firstChild of <name> element node.The are also
          two more text nodes, the '\ n' after <author> and </name>. So, <author> has 3 chids, tho text nodes '\ n'
          and <name> element node.
'''
import re, sys

class IloopsParser(object):
    '''
    Library used to parse iLoops xml outputs.
    Structure of the xml file (only the information considered):
    <protein>
       <name></name>
       <features>
           <feature>
               <loop_subClass></loop_subClass>
               <align_ini></align_ini>
               <align_end></align_end>
           </feature>
           [...]
       </features>
    </protein>
    [...]
    <interaction>
        <P1ID></P1ID>
        <P2ID></P2ID>
        <RF_prediction></RF_prediction>
        <RF_score></RF_score>
        <inferred_precisions>
            <inferred_precision>
                <unbalance_ratio>1</unbalance_ratio>
                <precision>0.74</precision>
                <error>0.015</error>
                </inferred_precision>
            <inferred_precision>
            [.x4.]
        </inferred_precisions>
        <interaction_signatures>
            <interaction_signature>
                <IS_type></IS_type>
                <P1IS></P1IS>
                <P2IS></P2IS>
                <pValue></pValue>
                </interaction_signature>
                [...]
        <interaction_signature>
    </interaction>
    [...]
    '''
    def __init__(self, xmlfile):
        '''
        Constructor that sets the xml file as a DOM object
        '''
        self.xmlfile = parse(xmlfile)

    def get_loops(self, protein_name, warnings=False):
        '''
        Get a list of the loops for a certain protein
        '''
        loops = []
        for protein in self.xmlfile.getElementsByTagName('protein'):
            if protein.getElementsByTagName('name')[0].firstChild.nodeValue == protein_name:
                for feature in protein.getElementsByTagName('feature'):
                    loop_subClass = feature.getElementsByTagName('loop_subClass')[0].firstChild.nodeValue
                    align_ini = int(feature.getElementsByTagName('align_ini')[0].firstChild.nodeValue)
                    align_end = int(feature.getElementsByTagName('align_end')[0].firstChild.nodeValue)
                    loops.append((loop_subClass, align_ini, align_end))
                continue
        if len(loops) == 0 and warnings:
            sys.stderr.write('No loops found for protein %s\n' % protein_name)
        return loops

    def get_interaction(self, protein1, protein2):
        '''
        Get a tuple of the interaction information (RF_prediction, RF_score, inferred_precision for ub 1, 
        inferred_precision for ub 10, inferred_precision for ub 20, inferred_precision for ub 50)
        '''
        for interaction in self.xmlfile.getElementsByTagName('interaction'):
            if interaction.getElementsByTagName('P1ID')[0].firstChild.nodeValue == protein1 and interaction.getElementsByTagName('P2ID')[0].firstChild.nodeValue == protein2:
                RF_prediction = interaction.getElementsByTagName('RF_prediction')[0].firstChild.nodeValue
                RF_score = float(interaction.getElementsByTagName('RF_score')[0].firstChild.nodeValue)
                inferred_precision_ubr1  = float(interaction.getElementsByTagName('inferred_precision')[0].childNodes[3].firstChild.nodeValue)
                inferred_precision_ubr10 = float(interaction.getElementsByTagName('inferred_precision')[1].childNodes[3].firstChild.nodeValue)
                inferred_precision_ubr20 = float(interaction.getElementsByTagName('inferred_precision')[2].childNodes[3].firstChild.nodeValue)
                inferred_precision_ubr50 = float(interaction.getElementsByTagName('inferred_precision')[3].childNodes[3].firstChild.nodeValue)
                return RF_prediction, RF_score, inferred_precision_ubr1, inferred_precision_ubr10, inferred_precision_ubr20, inferred_precision_ubr50
        sys.stderr.write('No interaction found for the protein tuple (%s,  %s)\n' % (protein1, protein2))
        return None

    def get_positive_signatures(self, protein1, protein2, max_protein_signatures = 3, warnings=False):
        '''
        Get the positive signatures of an interaction
        '''
        if max_protein_signatures > 3:
           max_protein_signatures = 3
        if max_protein_signatures < 1:
           max_protein_signatures = 1
        positive_signatures = []
        for interaction in self.xmlfile.getElementsByTagName('interaction'):
            if interaction.getElementsByTagName('P1ID')[0].firstChild.nodeValue == protein1 and interaction.getElementsByTagName('P2ID')[0].firstChild.nodeValue == protein2:
                for signature in interaction.getElementsByTagName('interaction_signature'):
                    if signature.childNodes[1].firstChild.nodeValue == 'positive':
                        P1IS = [x.strip('u\'') for x in re.split(', ', signature.childNodes[3].firstChild.nodeValue.strip('(,)'))]
                        P2IS = [x.strip('u\'') for x in re.split(', ', signature.childNodes[5].firstChild.nodeValue.strip('(,)'))]
                        pValue = float(signature.childNodes[7].firstChild.nodeValue)
                        if len(P1IS) <= max_protein_signatures and len(P2IS) <= max_protein_signatures:
                            positive_signatures.append((P1IS, P2IS, pValue))
        if len(positive_signatures) == 0 and warnings:
            sys.stderr.write('No positive signatures found for interaction (%s,  %s)\n' % (protein1, protein2))
        return positive_signatures

    def get_negative_signatures(self, protein1, protein2, max_protein_signatures = 3, warnings=False):
        '''
        Get the negative signatures of an interaction
        '''
        if max_protein_signatures > 3:
            max_protein_signatures = 3
        if max_protein_signatures < 1:
           max_protein_signatures = 1
        negative_signatures = []
        for interaction in self.xmlfile.getElementsByTagName('interaction'):
            if interaction.getElementsByTagName('P1ID')[0].firstChild.nodeValue == protein1 and interaction.getElementsByTagName('P2ID')[0].firstChild.nodeValue == protein2:
                for signature in interaction.getElementsByTagName('interaction_signature'):
                    if signature.childNodes[1].firstChild.nodeValue == 'negative':
                        P1IS = [x.strip('u\'') for x in re.split(', ', signature.childNodes[3].firstChild.nodeValue.strip('(,)'))]
                        P2IS = [x.strip('u\'') for x in re.split(', ', signature.childNodes[5].firstChild.nodeValue.strip('(,)'))]
                        pValue = float(signature.childNodes[7].firstChild.nodeValue)
                        if len(P1IS) <= max_protein_signatures and len(P2IS) <= max_protein_signatures:
                            negative_signatures.append((P1IS, P2IS, pValue))
        if len(negative_signatures) == 0 and warnings:
            sys.stderr.write('No negative signatures found for interaction (%s,  %s)\n' % (protein1, protein2))
        return negative_signatures

    def __str__(self):
        return self.xmlfile.toxml()

import puckepy
import unittest

#run in the current directory to run tests
#```
#$ python -m unittest test_all_modules.py
#```
#
# for verbose, run 
#```
#$ python -m unittest everything_test.py -v
#```
#
#
# Most tests are to see if the arguments from python match with their equivalent 
# rust function, see if all methods are available
#
# The setUp functions are tested implicitly, if anything fails here, that means there is an error 
# with the wrapping of the functions of the named parameter arguments.
# the self.setUp() method will not show in the methods tested, but will be called and its contents 
# are evaluated at runtime
#
#
#
#

# Test Formalism module
class Formalism(unittest.TestCase):

    def setUp(self) :

        self.fiveringPdb = puckepy.formalism.Pdb(filename="./data/fivering_adenosine.pdb").parse()
        self.sixringPdb = puckepy.formalism.Pdb(filename="./data/sixring_morpholino.pdb").parse()
        self.fiveringXyz = puckepy.formalism.Xyz(filename="./data/furanose.xyz").parse()

        # CP5() class
        self.cp51 = puckepy.formalism.CP5(0.5, 180.)
        self.cp52 = puckepy.formalism.CP5(amplitude=0.5, phase_angle=180.)
        self.cp53 = puckepy.formalism.CP5(0.5, 180.).invert()

        self.cp54 = puckepy.formalism.CP5().from_atomnames(
                pdb=self.fiveringPdb,
                query_names=["O4'", "C1'", "C2'", "C3'", "C4'", ]
                )

        self.cp55 = puckepy.formalism.CP5().from_indices(
                coordinates=self.fiveringPdb.coordinates,
                indices=[7, 8, 26, 24, 6]
                )

        # AS() class
        self.as1 = puckepy.formalism.AS(0.5, 180.)
        self.as2 = puckepy.formalism.AS(amplitude=0.5, phase_angle=180.)

        self.as3 = puckepy.formalism.AS().from_atomnames(
                pdb=self.fiveringPdb,
                query_names=["O4'", "C1'", "C2'", "C3'", "C4'", ]
                )

        self.as4 = puckepy.formalism.AS().from_indices(
                coordinates=self.fiveringPdb.coordinates,
                indices=[7, 8, 26, 24, 6]
                )


        # CP6() class
        self.cp61 = puckepy.formalism.CP6(0.5, 180., 90.)
        self.cp62 = puckepy.formalism.CP6(amplitude=0.5, phase_angle=180., theta=90.)
        self.cp63 = puckepy.formalism.CP6(0.5, 180., 90.).invert()

        self.cp64 = puckepy.formalism.CP6().from_atomnames(
                pdb=self.sixringPdb,
                query_names=["O5'", "C1'", "C2'", "N3'", "C4'", "C5'"]
                )

        self.cp65 = puckepy.formalism.CP6().from_indices(
                coordinates=self.sixringPdb.coordinates,
                indices=[6, 7, 27, 26, 23, 4]

                )

        # SP() class
        self.sp1 = puckepy.formalism.SP()

        self.sp2 = puckepy.formalism.SP().from_atomnames(
                pdb=self.sixringPdb,
                query_names=["O5'", "C1'", "C2'", "N3'", "C4'", "C5'"]
                )

        self.sp3 = puckepy.formalism.SP().from_indices(
                coordinates=self.sixringPdb.coordinates,
                indices=[6, 7, 27, 26, 23, 4]
                )

    def test_formalism(self):
        self.assertEqual('foo'.lower(), 'foo')


    # CP5 Testings
    @unittest.expectedFailure
    def test_incorrect_amount_of_indices_CP5(self):
        self.cp5a = puckepy.formalism.CP5().from_indices(
                coordinates=self.fiveringPdb.coordinates,
                indices=[7, 8, 26, 24, 6, 2]
                )

    @unittest.expectedFailure
    def test_incorrect_amount_of_querynames_CP5(self):
        self.cp5b = puckepy.formalism.CP5().from_atomnames(
                pdb=self.fiveringPdb,
                query_names=["O4'", "C1'", "C2'", "C3'", "C4'", "C5'"]
                )

    # AS Testings
    @unittest.expectedFailure
    def test_incorrect_amount_of_indices_AS(self):
        self.asa = puckepy.formalism.AS().from_indices(
                coordinates=self.fiveringPdb.coordinates,
                indices=[7, 8, 26, 24, 6, 2]
                )

    @unittest.expectedFailure
    def test_incorrect_amount_of_querynames_AS(self):
        self.asb = puckepy.formalism.AS().from_atomnames(
                pdb=self.fiveringPdb,
                query_names=["O4'", "C1'", "C2'", "C3'", "C4'", "C5'"]
                )

    # CP6 Testings
    @unittest.expectedFailure
    def test_incorrect_amount_of_indices_CP6(self):
        self.cp6a = puckepy.formalism.CP6().from_indices(
                coordinates=self.sixringPdb.coordinates,
                indices=[7, 8, 26, 24, 6]
                )

    @unittest.expectedFailure
    def test_incorrect_amount_of_querynames_CP6(self):
        self.cp6b = puckepy.formalism.CP6().from_atomnames(
                pdb=self.sixringPdb,
                query_names=["O4'", "C1'", "C2'", "C3'", "C4'"]
                )

    # SP Testings
    @unittest.expectedFailure
    def test_incorrect_amount_of_indices_SP(self):
        self.sp6a = puckepy.formalism.SP().from_indices(
                coordinates=self.sixringPdb.coordinates,
                indices=[7, 8, 26, 24, 6]
                )

    @unittest.expectedFailure
    def test_incorrect_amount_of_querynames_SP(self):
        self.sp6b = puckepy.formalism.SP().from_atomnames(
                pdb=self.sixringPdb,
                query_names=["O4'", "C1'", "C2'", "C3'", "C4'"]
                )
# Test Confsampling module
class Confsampling(unittest.TestCase):

    def setUp(self) :
        self.a1 = puckepy.confsampling.Sixring(amount=630)
        self.a2 = puckepy.confsampling.SixringAxes(amount=630)

        self.b1 = puckepy.confsampling.Fivering(interval=21)
        self.b2 = puckepy.confsampling.FiveringAxes(interval=21)

        self.c1 = puckepy.confsampling.Peptide(interval=37)
        self.c2 = puckepy.confsampling.PeptideAxes(interval=37)

    
    def test_confsampling(self):
        self.assertEqual('foo'.lower(), 'foo')

# Test Geometry module
class Geometry(unittest.TestCase):


    def setUp(self) :
        self.fiveringXyz = puckepy.formalism.Xyz(filename="./data/furanose.xyz").parse()


        self.dihedral = puckepy.geometry.dihedral(
                p0= self.fiveringXyz[0],
                p1= self.fiveringXyz[1],
                p2= self.fiveringXyz[2],
                p3= self.fiveringXyz[3],
                )
        self.bondangle = puckepy.geometry.bondangle(
                p0= self.fiveringXyz[0],
                p1= self.fiveringXyz[1],
                p2= self.fiveringXyz[2],
                )
        self.bondlength = puckepy.geometry.bondlength(
                p0= self.fiveringXyz[0],
                p1= self.fiveringXyz[1],
                )

    def test_geometry(self):
        self.assertEqual('foo'.lower(), 'foo')

if __name__ == '__main__':
    unittest.main()

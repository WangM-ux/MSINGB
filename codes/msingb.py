import argparse
import sys
import os

from getModel import model_pre
from getTag import Test_Input

datapath = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
dataPath = os.path.abspath(os.path.dirname(os.getcwd()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage="it's usage tip.",  description="Accept MAF file in response to user")
    parser.add_argument("--maf", required=True, help="input maf file")
    parser.add_argument("--outdir", required=True, help="input out dir")
    args = parser.parse_args()
    input = args.maf
    output = args.outdir

    srsPath = dataPath + '/data/simpleRepeat.txt'
    trainPath = dataPath + '/data/all_train.csv'
    tf = Test_Input(inputMaf=input, output=output)
    feature = tf.creat_tagged_feature(ref_repeats_file=srsPath)
    prediction = model_pre(trainSet=trainPath, inputFeature=feature, output=output)




import argparse
import sys
import os
curPath = os.path.abspath(os.path.dirname(__file__))
rootPath = os.path.split(curPath)[0]
sys.path.append(rootPath)
from MSINGB_pre.getModel import model_pre
from MSINGB_pre.getTag import Test_Input

if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage="it's usage tip.",  description="Accept MAF file in response to user")
    parser.add_argument("--maf", required=True, help="input maf file")
    parser.add_argument("--outdir", required=True, help="input out dir")
    args = parser.parse_args()
    input = args.maf
    output = args.outdir

    srsPath = curPath + '../data/simpleRepeat.txt'
    trainPath = curPath + '../data/all_train.csv'
    tf = Test_Input(inputMaf=input, output=output)
    feature = tf.creat_tagged_feature(ref_repeats_file=srsPath)
    prediction = model_pre(trainSet=trainPath, inputFeature=feature, output=output)




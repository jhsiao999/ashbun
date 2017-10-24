Many methods for single cell abundance analysis
=====================================================




# Structure of functions
for a single count table, run
|---- query.evaluation(): generate results of false discovery rate and true positive rate
      |-- query.pipeline(): combine results of multiple DE and multiple normalization methods
          |-- filter.Wrapper(): filter genes and samples
          |-- query.methodsMeanExpression(): run multiple DE methods
              |-- methodWrapper.DESeq2()
              |-- methodWrapper.edgeR()
              |-- methodWrapper.limmaVoom()
              |-- methodWrapper.bpsc
              |-- methodWrapper.rots
              |-- methodWrapper.scde
          |-- query.methodsNormalization(): run multiple normalization methods
              |-- methodsNormalize.R

function j= modelJacobian_weberWrapper(t,y,N0,growthOn,diffusionOn)
  
  [modelParams,growth]= singlecell.modelCoeffs_weber(y,N0,growthOn,diffusionOn);
  j= singlecell.modelJacobian_weber(t,y,modelParams,growth);
end

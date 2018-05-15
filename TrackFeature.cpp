void R2Image::TrackFeatures(R2Image *otherImage){
  


  vector<vector<double> > matching_features;

  
  int searchAreaHalfWidth = otherImage->width * 0.01;
  int searchAreaHalfHeight = otherImage->height * 0.01;
   cout << searchAreaHalfWidth << endl;
   cout << searchAreaHalfHeight << endl;

  //search for all feature in vector
  for (int feature = 0; feature < imgA_features.size(); feature++) {
    double minSum = 100000; 
    double minX = 0;
    double minY = 0;
    for (int i = -searchAreaHalfWidth; i <=searchAreaHalfWidth;i++) {
      for (int j = -searchAreaHalfHeight; j <= searchAreaHalfHeight; j++) {

        if (imgA_features[feature][1]+i >= 0 && imgA_features[feature][1]+i< width &&
                imgA_features[feature][2]+j>=0 && imgA_features[feature][2]+j<height){

                 // start a window
                  R2Pixel *ssd = new R2Pixel();
                  for (int k = -3; k <=3; k++) {
                    for (int l = -3; l <=3; l++) {
                      if (imgA_features[feature][1]+k >= 0 && imgA_features[feature][1]+k < otherImage->width &&
                          imgA_features[feature][2]+l >=0 && imgA_features[feature][2]+l < otherImage->height &&
                          imgA_features[feature][1]+i+k >= 0 && imgA_features[feature][1]+i+k < otherImage->width &&
                          imgA_features[feature][2]+j+l>=0 && imgA_features[feature][2]+j+l<otherImage->height) 
             
                          *ssd += (Pixel(imgA_features[feature][1]+k,imgA_features[feature][2]+l)-otherImage->Pixel(imgA_features[feature][1]+i+k,imgA_features[feature][2]+j+l))*(Pixel(imgA_features[feature][1]+k,imgA_features[feature][2]+l)-otherImage->Pixel(imgA_features[feature][1]+i+k,imgA_features[feature][2]+j+l));          
                      
                    }
                  }
                  double sum = ssd->Red()+ssd->Green()+ssd->Blue();

                  if (sum < minSum) {
                    minSum = sum; 
                    minX = imgA_features[feature][1]+i;
                    minY = imgA_features[feature][2]+j;
                  
                  }

        }
            
      }
    }
       cout << "difference: "<< minSum << endl;
    
      vector<double> item;
      item.push_back(minX);
      item.push_back(minY);
      matching_features.push_back(item);
    
  }

  //draw the matching features 
  for (int i = 0; i< matching_features.size(); i++){
    line(imgA_features[i][1],matching_features[i][0],imgA_features[i][2],matching_features[i][1],0,1,0);
  }

  imgA_features.clear();

  for(int i = 0; i < matching_features.size();i++){
    vector<double> object;
    object.push_back(0.0);
    // cout << "matching features: "<< matching_features[i][0]<< " "<<matching_features[i][1]<<endl;
    object.push_back(matching_features[i][0]);
    object.push_back(matching_features[i][1]);

    imgA_features.push_back(object);
  }

   



}
String[] lines, lines2, lines3, lines4,lines5, pieces, pieces2, pieces3, pieces4,pieces5;
int i,j,jmax,k;
int ifac = 1;
float col1= 0;
float col2= 0;
float Aux_ave=0;
float Cyt_ave=0;
int cell_type=0;
int cell_lineage=0;
int cell_defect=0;

//float a=8,b=220,c=210,d=2;//initial3,4,5,8,9,10
float a=8,b=220,c=240,d=2;//initial1,2,6

float g=10000000;
float af=1000000,bf=100;
int arrow_headlen=10;
int arrow_angle=10;
FloatList x = new FloatList ();
FloatList y = new FloatList ();
FloatList Fx = new FloatList ();
FloatList Fy = new FloatList ();
float x_ch=0;
float y_ch=0;
float S1=0;
float S2=0;
float Theta1=0;
float Theta2=0;

float Stress_ave=0;
float Stress_vari=0;
float Anisotoropy_ave=0;
float Va_cellnumber=0;
void setup() {
  //size(750,750);
  size(850,850);
  background(255);
  strokeWeight(3);
  frameRate(100);//100
  smooth();
}

void draw() {
  if(ifac<2532/*2782*//*2532*//*3282*/) {
   background(255);
   strokeWeight(3);
   //vertex
   //lines = loadStrings("/Users/fujiwara/Desktop/vascul_model/data/z_cellform"+ifac+".dat");
   //A-C
   lines  = loadStrings("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/data_form2/z_cellform"+ifac+".dat");
   lines3 = loadStrings("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/data_form2/z_celltype"+ifac+".dat");
   lines4 = loadStrings("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/data_form2/z_stresstensor"+ifac+".dat");
   lines5 = loadStrings("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/data_form2/z_centroid"+ifac+".dat");
   
   /*lines  = loadStrings("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/sketch_vasbund77/initial2/hanx6_ini2_707/data_form2/z_cellform"+ifac+".dat");
   lines3 = loadStrings("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/sketch_vasbund77/initial2/hanx6_ini2_707/data_form2/z_celltype"+ifac+".dat");
   lines4 = loadStrings("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/sketch_vasbund77/initial2/hanx6_ini2_707/data_form2/z_stresstensor"+ifac+".dat");
   lines5 = loadStrings("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/sketch_vasbund77/initial2/hanx6_ini2_707/data_form2/z_centroid"+ifac+".dat");
   */
   //pieces3 = split(lines2[0],  '\t');
  // Aux_ave = float(pieces3[0]);
  // Cyt_ave = float(pieces3[1]);
   
/*for(i=0; i<lines4.length; i++)//initial_va
{
  if (cell_type==1 || cell_type==2 || cell_type==3 || cell_type==6) {
     pieces4 = split(lines4[i],  ' ');
     S1=float(pieces4[0]);
     S2=float(pieces4[1]);
    Stress_ave+=S1+S2;
     //Anisotoropy_ave+=S1-S2;
     Va_cellnumber+=1;
  }
}
//Stress_ave=Stress_ave/lines4.length;
  Stress_ave=Stress_ave/Va_cellnumber;
//Anisotoropy_ave=Anisotoropy_ave/lines4.length;*/

/*Stress_ave=0;
for(i=18; i<26; i++)//initial_va endodermis
{
     pieces4 = split(lines4[i],  ' ');
     S1=float(pieces4[0]);
     S2=float(pieces4[1]);
    Stress_ave+=S1+S2;
}
Stress_ave=Stress_ave/8;

Stress_vari=0;
for(i=18; i<26; i++)//initial_va
{
     pieces4 = split(lines4[i],  ' ');
     S1=float(pieces4[0]);
     S2=float(pieces4[1]);
    Stress_vari+=(S1+S2-Stress_ave)*(S1+S2-Stress_ave);
}
Stress_vari=sqrt(Stress_vari/8);

 Anisotoropy_ave=0;
 Va_cellnumber=0;
 for(i=0; i<lines.length; i++)//initial_va
  {
   pieces3 = split(lines3[i],  '\t');
   cell_type=int(pieces3[0]);
 if (cell_type==1 || cell_type==2 || cell_type==3 || cell_type==6) {
    pieces4 = split(lines4[i],  ' ');
     S1=float(pieces4[0]);
     S2=float(pieces4[1]);
   Anisotoropy_ave+=S1-S2;
       Va_cellnumber+=1;
  // }
  }
 Anisotoropy_ave=Anisotoropy_ave/Va_cellnumber;*/


  for(i=0; i<lines.length; i++)//initial_va
  {
     pieces = split(lines[i],  '\t');
    for (j=0; j<(pieces.length-1)*0.5; j++)
    {
    x.append(b+(float(pieces[2*j]))*a) ;
    y.append(c+(float(pieces[2*j+1]))*a);
    }
    
   pieces3 = split(lines3[i],  '\t');
   cell_type=int(pieces3[0]);
   cell_lineage=int(pieces3[1]);
   cell_defect=int(pieces3[2]);
   // pieces2 = split(lines2[i+1],  '\t'); 
   // col1=(float(pieces2[0])/Aux_ave);
   // col2=(float(pieces2[1])/Cyt_ave);
   // println(col1,col2);
   // stroke(255);
   // fill(0,0,255,200);
   
    pieces4 = split(lines4[i],  ' ');
     S1=float(pieces4[0]);
     S2=float(pieces4[1]);
     Theta1=float(pieces4[2]);
     Theta2=float(pieces4[3]);
     
    stroke(0);
    if (cell_type==1) {
    fill(255,255,0);//yellow-xylem
    // fill(139,69,19);
   // noFill();
   }
    if (cell_type==2) {
    fill(12,0,200);//blue-pse
     //fill(255,0,0);
    // noFill();
   }
   if (cell_type==7) {
    fill(12,0,200);//blue-pse-ln
     //fill(255,0,0);
   //  noFill();
   }
   
    if (cell_type==3) {
    fill(12,150,220);//aqua-ipc
     //fill(0,255,0);
   //  noFill();
   }
    if (cell_type==4) {
    fill(128,128,128);//gray-peri
   }
    if (cell_type==5) {
    fill(202,204,178);//lightgold-end
    //fill(190,162,202);
   }
    if (cell_type==6) {
   fill(150,12,220);//purple-opc
    //fill(255,255,0);
   // noFill();
   }
   
   /*if (cell_type==1) {
    fill(255,255,0);
   }
    if (cell_type==2) {
    fill(12,0,200);
   }
    if (cell_type==3) {
    fill(25,135,220);
   }
    if (cell_type==4) {
    fill(140,140,140);
   }
    if (cell_type==5) {
    fill(202,204,178);
    //fill(190,162,202);
   }
    if (cell_type==6) {
    //fill(12,100,200);
    fill(25,135,220);
   }
    if (cell_type==7) {
    //fill(12,100,200);
    fill(200,135,220);
   }*/
   /*if(cell_lineage==5){
      //stroke(0,100);
      fill(200,0,0);
      //noFill();
      //noStroke();
    }
    if(cell_lineage==6){
      //stroke(0,100);
      fill(0,200,0);
      //noFill();
      //noStroke();
    }
    if(cell_lineage==7){
      //stroke(0,100);
      fill(0,200,200);
      //noFill();
      //noStroke();
    }
    if(cell_lineage==8){
      //stroke(0,100);
      fill(200,100,100);
      //noFill();
      //noStroke();
    }*/
    
    //OPC defect
   if(cell_defect==1){
      //stroke(0,100);
      fill(100,0,0);
      //noFill();
      //noStroke();
    }
    if(cell_lineage==7){
    //if(i==14){
      //stroke(0,100);
      //fill(200,100,100);
      //noFill();
      //noStroke();
    }
    /* if (cell_type==6) {
   fill(150,12,220);//purple
    //fill(255,255,0);
   }*/
    
    //noFill();
    //fill(5*cell_lineage,0,0);
   
   //stress 
   if (cell_type==1||cell_type==2||cell_type==3||cell_type==6||cell_type==7) {
     //fill(0.5*(S1+S2)*g);
      }
     // println(S1+S2);
     //fill(50*(S1+S2)/(Stress_ave)+50);
     
     // fill(int(10*sqrt((S1+S2-Stress_ave)*(S1+S2-Stress_ave))/(Stress_vari)+50));
      //fill(8000*(S1+S2)+50);
      // println(sqrt(pow(S1+S2-Stress_ave,2))/(Stress_vari));
    // println(int(10*sqrt((S1+S2-Stress_ave)*(S1+S2-Stress_ave))/(Stress_vari)+50));
    
 //noFill();
 strokeWeight(3);
  beginShape();
  for (j=0; j<x.size(); j++)
  {
  vertex((float)x.get(j),(float)y.get(j));
  //println((float)x.get(j),(float)y.get(j));
  }
  vertex((float)x.get(0),(float)y.get(0));
  endShape();
  jmax=x.size();
   for (j=0; j<jmax; j++)
  {
  x.remove(0);
  y.remove(0);
  }
  
  //centroid
     pieces5 = split(lines5[i],  ' ');
     x_ch=(b+(float(pieces5[0]))*a);
     y_ch=(c+(float(pieces5[1]))*a);
      stroke(220);
      strokeWeight(4);
      if (cell_type==1 || cell_type==2 || cell_type==3 || cell_type==7 || cell_type==6) {
      if (S1-S2>Anisotoropy_ave*0.3) {
    //line(x_ch-(S1-S2)*cos(Theta1)*0.5/Anisotoropy_ave*30,y_ch-(S1-S2)*sin(Theta1)*0.5/Anisotoropy_ave*30,x_ch+(S1-S2)*cos(Theta1)*0.5/Anisotoropy_ave*30,y_ch+(S1-S2)*sin(Theta1)*0.5/Anisotoropy_ave*30);
   //line(x_ch-(S1-S2)*cos(Theta1)*0.5*g,y_ch-(S1-S2)*sin(Theta1)*0.5*g,x_ch+(S1-S2)*cos(Theta1)*0.5*30,y_ch+(S1-S2)*sin(Theta1)*0.5*g);  
    }
      }
      
     if (ifac==2781/*2781*//*2531*//*3281*/) {
        if (cell_type==1) {
       //println(float(pieces5[0]),float(pieces5[1]),S1,S2,Theta1,Theta2);
        }
      }
      if (ifac==2781/*2781*//*2531*//*3281*/) {
        if (cell_type==3 || cell_type==6) {
       //println(float(pieces5[0]),float(pieces5[1]),S1,S2,Theta1,Theta2);
        }
      }
      if (ifac==2531/*2781*//*2531*//*3281*/) {
        if (cell_type==2) {
       //println(float(pieces5[0]),float(pieces5[1]),S1,S2,Theta1,Theta2);
       println(float(pieces5[0]),float(pieces5[1]));
        }
      }
      
   fill(0);
  textSize(15);
     //number_character(cell_type,x_ch-10,y_ch);
   //number_character(i,x_ch-10,y_ch);
  } //i_cell
  
   for(i=0; i<lines.length; i++) {//pse_edge line
   pieces3 = split(lines3[i],  '\t');
   cell_type=int(pieces3[0]);
    if (cell_type==2) {
     pieces = split(lines[i],  '\t');
    for (j=0; j<(pieces.length-1)*0.5; j++) {
    x.append(b+(float(pieces[2*j]))*a) ;
    y.append(c+(float(pieces[2*j+1]))*a);
    }
   stroke(220,0,0);
   noFill();
  // noStroke();
  beginShape();
  for (j=0; j<x.size(); j++) {
  vertex((float)x.get(j),(float)y.get(j));
  //println((float)x.get(j),(float)y.get(j));
  }
  vertex((float)x.get(0),(float)y.get(0));
  endShape();
  jmax=x.size();
   for (j=0; j<jmax; j++) {
    x.remove(0);
    y.remove(0);
    }
    }
   }
 
/*
for(i=0; i<lines.length; i++) {//xylem_edge line
   pieces3 = split(lines3[i],  '\t');
   cell_type=int(pieces3[0]);
    if (cell_type==1) {
     pieces = split(lines[i],  '\t');
    for (j=0; j<(pieces.length-1)*0.5; j++) {
    x.append(b+(float(pieces[2*j]))*a) ;
    y.append(c+(float(pieces[2*j+1]))*a);
    }
   stroke(255,255,0);
   noFill();
  // noStroke();
  beginShape();
  for (j=0; j<x.size(); j++) {
  vertex((float)x.get(j),(float)y.get(j));
  //println((float)x.get(j),(float)y.get(j));
  }
  vertex((float)x.get(0),(float)y.get(0));
  endShape();
  jmax=x.size();
   for (j=0; j<jmax; j++) {
    x.remove(0);
    y.remove(0);
    }
    }
   }*/

/*for(i=27; i<lines.length; i++)//divi
  {
     pieces = split(lines[i],  '\t');
    for (j=0; j<(pieces.length-1)*0.5; j++)
    {
    x.append(b+(float(pieces[2*j]))*a) ;
    y.append(c+(float(pieces[2*j+1]))*a);
    }
    
   // pieces2 = split(lines2[i+1],  '\t'); 
   // col1=(float(pieces2[0])/Aux_ave);
   // col2=(float(pieces2[1])/Cyt_ave);
   // println(col1,col2);
   // stroke(0);
   // fill(int(125*col2),int(125*col1),0);
    //stroke(255);
    //fill(200,100,255,200);
    stroke(255,0,0);
   noFill();
  beginShape();
  for (j=0; j<x.size(); j++)
  {
  vertex((float)x.get(j),(float)y.get(j));
  //println((float)x.get(j),(float)y.get(j));
  }
  vertex((float)x.get(0),(float)y.get(0));
  endShape();
  jmax=x.size();
   for (j=0; j<jmax; j++)
  {
  x.remove(0);
  y.remove(0);
  }
  }*/
  
 /* lines4 = loadStrings("/Users/fujiwara/Desktop/cellform_c++/za2_force"+ifac+".dat");
  for(i=0; i<lines4.length; i++)
  {
     pieces4 = split(lines4[i],  '\t');
    for (j=0; j<pieces4.length; j++)
    {
    x.append(b+(float(pieces4[0]))*a) ;
    y.append(c+(float(pieces4[1]))*a);
    Fx.append((float(pieces4[2]))*af) ;
    Fy.append((float(pieces4[3]))*af);
    }
    stroke(255,0,0);
 
  for (j=0; j<x.size(); j++)
  {
  arrowline((float)x.get(j), (float)y.get(j), (float)x.get(j)+(float)Fx.get(j), (float)y.get(j)+(float)Fy.get(j));
  }
  jmax=x.size();
   for (j=0; j<jmax; j++)
  {
  x.remove(0);
  y.remove(0);
  Fx.remove(0);
  Fy.remove(0);
  }
  }
  */
  println(ifac);
  }
  
  //fill(0);
 //text("(c)2017 Motohiro Fujiwara", 540,650);
 
/* for (i=0;i<10;i++){
 fill(25*i);
 noStroke();
 rect(100+20*i,100,20,20);
 }*/
 
/* strokeWeight(1);
 for( j = 0; j < 100; j++){
  stroke(255 * j*0.01);
  line(j+100, 100, j+100, 120);
}
  stroke(0);
  line(100, 95, 100, 100);
  line(150, 95, 150, 100);
  line(200, 95, 200, 100);*/

   ifac+=10;//5
   //saveFrame("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/sketch_vasbund77/initial10/partiallyaniso_para/partaniso/z_form_####.png");
 
   if (ifac==61 /*2781*//*2531*//*3281*/) {
    //saveFrame("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/sketch_vasbund77/initial2/z_form_####.png");
    //saveFrame("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/sketch_vasbund77/z_form_1.png");
     }
   if (ifac==3281/*2781*//*2531*//*3281*/) {
    //saveFrame("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/sketch_vasbund77/initial2/z_form_####.png");
    //saveFrame("/Users/fujiwara/Desktop/vertexmodel/vascularsimu/sketch_vasbund77/z_form_3281.png");
     }
     
   
/*if (ifac > 3400) {
lines4 = loadStrings("/Users/fujiwara/Desktop/vascul_model/result2/cell-centroid_result.dat");
  for(i=0; i<lines4.length; i++)  {
     pieces4 = split(lines4[i],  ' ');
     x_ch=(b+(float(pieces4[0]))*a);
     y_ch=(c+(float(pieces4[1]))*a);
     fill(0);
     textSize(15);
     number_character(i,x_ch-10,y_ch);
  }
}*/

}

void arrowline(float x1, float y1, float x2, float y2){  
  line(x1, y1, x2, y2);
  float dx = x2-x1;
  float dy = y2-y1;
  float theta = atan2(dy, dx)-PI;
  arrow_headlen = (int)(0.3*sqrt(dx*dx+dy*dy));
  line(x2, y2, arrow_headlen*cos(theta+radians(arrow_angle))+x2, arrow_headlen*sin(theta+radians(arrow_angle))+y2);
  line(x2, y2, arrow_headlen*cos(theta-radians(arrow_angle))+x2, arrow_headlen*sin(theta-radians(arrow_angle))+y2);
   
}
void number_character(int i,float x, float y){
 text(i, x,y);
}

void mousePressed() {
  loop();  // Holding down the mouse activates looping
}

void mouseReleased() {
  noLoop();  // Releasing the mouse stops looping draw()
}

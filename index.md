# Introduction to Bulk RNA-seq data analysis
### 17th, 20th and 21st March 2023
#### Taught Hybrid
#### Bioinformatics Training Facility, Craik-Marshall Building, Downing Site, University of Cambridge

![](Bulk_RNAseq_Course_Base/images/CRUK_Cambridge_Major_Centre_logo.jpg)

## Instructors

* Abigail Edwards - Bioinformatics Core, Cancer Research UK Cambridge Institute
* Ashley D Sawle - Bioinformatics Core, Cancer Research UK Cambridge Institute
* Sina Beier - MRC Toxicology Unit, University of Cambridge
* Raquel Manzano Garcia - Cancer Research UK Cambridge Institute
* Jonathon Price - Department of Biochemistry, University of Cambridge
* Jiayin Hong - Department of Biochemistry, University of Cambridge

**Helpers**
* Ulrika Yuan - Department of Biochemistry, University of Cambridge


## Outline

In this workshop, you will be learning how to analyse RNA-seq data. This will
include read alignment, quality control, quantification against a reference,
reading the count data into R, performing differential expression analysis, and
gene set testing, with a focus on the DESeq2 analysis workflow. You will learn
how to generate common plots for analysis and visualisation of gene expression
data, such as boxplots and heatmaps.

This workshop is aimed at biologists interested in learning how to perform
differential expression analysis of RNA-seq data.

Whilst we have run this course for several years, we are still improving and learning how to
teach it best.  
Please bear with us if there are any technical hitches, and
be aware that timings for different sections laid out in the schedule below may
not be adhered to. There may be some necessity to make adjustments to the course
as we go.

> ## Prerequisites
>
> __**Some basic experience of using a UNIX/LINUX command line is assumed**__
>
> __**Some R knowledge is assumed and essential. Without it, you
> will struggle on this course.**__
> If you are not familiar with the R statistical programming language we
> strongly encourage you to work through an introductory R course before
> attempting these materials.
> We recommend our [Introduction to R course](https://bioinformatics-core-shared-training.github.io/r-intro/)

## Shared Google Document

This
[Google Document](https://docs.google.com/document/d/1zClM_fTWQIECQUkzWSBGm66s0Giv6qaSkfM3KyCVd0c/edit?usp=sharing) contains useful information and links..

Please use it to post any long form questions you have during the course.

The trainers will be monitoring the document and will answer questions as quickly
as they can.

## Introduce Yourself

There is another Google Doc
[Google Document](https://docs.google.com/document/d/1GGAJlyF45s8MyZjSIWh8JeC9AJ8DBA131D8kYRsM0h8/edit?usp=sharing)
Please write a couple sentences here to introduce yourself to the class, tell
us a bit about your background and what you hope to get out of this course.  If
you are a student or staff at the University of Cambridge, tell us which
Department you are in.


## Course etiquette

As this course is being partially taught online and there are a large number of participants,
online students will all need to follow a [few simple rules](https://docs.google.com/presentation/d/e/2PACX-1vQv9nTlsdRC9iZJU138tLL1jrwNoryp8P-FnXxb_ugOOWjbav4QHTLYLLZj2KK4kTO0_3x3VlzSdrUu/pub?start=false&loop=false&delayms=3000) to ensure things run as smoothly as possible:

1. Please mute your microphone

2. To get help from a tutor, please click the "Raise Hand" button in Zoom:

    ![](Bulk_RNAseq_Course_Base/images/raise_hand.png)

   This can be found by clicking on the 'Participants' button. A tutor will
   then contact you in the chat. If necessary, you and the tutor can be moved
   to a breakout room where you can discuss your issue in more detail.

3. Please ask any general question by typing it into the Google Doc mentioned above

4. During practicals, when you are done, please press the green "Yes" button:

    ![](Bulk_RNAseq_Course_Base/images/yes_button.png)

   This way we will know when we can move on.
   
trainers will be monitoring the zoom at all times and communicating with the presenters onsite and will be able to pass on any messages/questions.

## Timetable

**This is the first time we have offered this course as a hybrid class, all times here should be
regarded as aspirations**

### Day 1

**Trainers:** Abbi, Ash, Jiayin, Sina, Jon (OL), Raquel (OL), Ulrika (OL)

9:30 - 9:45 - Welcome! <!-- Abbi -->

9:45 - 10:15 - [Introduction to RNAseq
Methods](Bulk_RNAseq_Course_Base/Markdowns/01_Introduction_to_RNAseq_Methods.html)- Sina Beier

10:15 - 11:00 [Raw read file format and
QC](Bulk_RNAseq_Course_Base/Markdowns/02_FastQC_introduction.html) - Sina Beier
    - [Practical](Bulk_RNAseq_Course_Base/Markdowns/02_FastQC_practical.html)
<!--    - [Practical solutions](Bulk_RNAseq_Course_Base/Markdowns/02_FastQC_practical.Solutions.html)  -->

11:00 - 13:30 [Alignment and Quantification of Gene Expression with Salmon](Bulk_RNAseq_Course_Base/Markdowns/03_Quantification_with_Salmon_introduction.html) - Ashley Sawle  
    - [Practical](Bulk_RNAseq_Course_Base/Markdowns/03_Quantification_with_Salmon_practical.html)  
    - [Practical solutions](Bulk_RNAseq_Course_Base/Markdowns/03_Quantification_with_Salmon_practical.Solutions.html)
    
13:30 - 14:30 Lunch

14:30 - 15:30 [QC of alignment](Bulk_RNAseq_Course_Base/Markdowns/04_Quality_Control_introduction.html) - Abbi Edwards

  - [Practical](Bulk_RNAseq_Course_Base/Markdowns/04_Quality_Control_practical.html) ([pdf](Bulk_RNAseq_Course_Base/Markdowns/04_Quality_Control_practical.pdf))  
<!--    - [Practical solutions](Bulk_RNAseq_Course_Base/Markdowns/04_Quality_COntrol_practical.Solutions.html) ([pdf](Bulk_RNAseq_Course_Base/Markdowns/04_Quality_Control_practical.Solutions.pdf))    -->

15.30 - 17.30 [Data Exploration in R](Bulk_RNAseq_Course_Base/Markdowns/05_Data_Exploration.html) ([pdf](Bulk_RNAseq_Course_Base/Markdowns/05_Data_Exploration.pdf)) - Jiayin Hong   
- [Practical solutions](Bulk_RNAseq_Course_Base/Markdowns/05_Data_Exploration.Solutions.html) ([pdf](Bulk_RNAseq_Course_Base/Markdowns/05_Data_Exploration.Solutions.pdf))  
- [Live script](Live_Scripts/data_exploration.R) 

### Day 2

**Trainers:** Abbi, Ash, Jiayin, Jon, Sina (OL), Raquel (OL), Ulrika (OL)

<!-- Welcome Announcements - Abbi -->
9:30 - 10:15  [Introduction to RNAseq Analysis in
R](Bulk_RNAseq_Course_Base/Markdowns/06_Introduction_to_RNAseq_Analysis_in_R.html) - Jiayin Hong

10:15 - 13:00 Statistical Analysis of Bulk RNAseq Data

- Part I: [Statistics of RNA-seq analysis](Bulk_RNAseq_Course_Base/additional_scripts_and_materials/RNA-seq_stats.pdf) - Abbi Edwards

- Part II: [Linear Models in R and DESeq2](Bulk_RNAseq_Course_Base/Markdowns/07_Linear_Models.html) ([pdf](Bulk_RNAseq_Course_Base/Markdowns/07_Linear_Models.pdf)) - Abbi Edwards  
    - [Slides](Bulk_RNAseq_Course_Base/additional_scripts_and_materials/Statistical_models_in_R_DESeq2.pdf) 
    - Find the worksheet in `Course_Materials/stats/models_in_r_worksheet.R`

13:00 - 14:00 Lunch

14:00 - 17:30 - [Differential Expression for RNA-seq](Bulk_RNAseq_Course_Base/Markdowns/08_DE_analysis_with_DESeq2.html)
([pdf](Bulk_RNAseq_Course_Base/Markdowns/08_DE_analysis_with_DESeq2.pdf)) - Jon Price
   - [practical solutions](Bulk_RNAseq_Course_Base/Markdowns/08_DE_analysis_with_DESeq2.Solutions.html) ([pdf](Bulk_RNAseq_Course_Base/Markdowns/08_DE_analysis_with_DESeq2.Solutions.pdf)) 
  - [live script](Live_Scripts/DESEQ2SESSION.R)
<!--      - [extra models plots](Bulk_RNAseq_Course_Base/additional_scripts_and_materials/Expl_all.pdf) -->

### Day 3

**Trainers:** Abbi, Ash, Raquel, Sina (OL), Jon (OL), Jiayin (OL)

9.30 - 9.45 - [Recap of Day 1 and 2](Bulk_RNAseq_Course_Base/additional_scripts_and_materials/Analysis_of_RNA-seq_data_day3recap.pdf) - Abbi Edwards

9.45 - 12.30 [Annotation and Visualisation of RNA-seq results](Bulk_RNAseq_Course_Base/Markdowns/09_Annotation_and_Visualisation.html) - Raquel Manzano Garcia
 <!-- - [practical solutions](Bulk_RNAseq_Course_Base/Markdowns/09_Annotation_and_Visualisation_Solutions.html)   -->
 <!-- - [live script](live_scripts/AandV.R) -->

12.30 - 13.30 Lunch

13.30 - 16:30  [Gene-set testing](Bulk_RNAseq_Course_Base/Markdowns/10_Gene_set_testing_introduction.html) - Abbi Edwards
   - [Practical (html)](Bulk_RNAseq_Course_Base/Markdowns/10_Gene_set_testing.html) [(rmd)](Bulk_RNAseq_Course_Base/Markdowns/10_Gene_set_testing.Rmd) [(pdf)](Bulk_RNAseq_Course_Base/Markdowns/10_Gene_set_testing.pdf)
  <!-- - [Practical solutions (html)](Bulk_RNAseq_Course_Base/Markdowns/10_Gene_set_testing.Solutions.html) [(rmd)](Bulk_RNAseq_Course_Base/Markdowns/10_Gene_set_testing.Solutions.Rmd) [(pdf)](Bulk_RNAseq_Course_Base/Markdowns/10_Gene_set_testing.Solutions.pdf) -->
  <!-- - [Live Script](live_scripts/GS.R) -->

<!-- Goodbye: Abbi -->

## Source Materials for Practicals

The lecture slides and other source materials, including R code and
practical solutions, can be found in the course's [Github
repository](https://github.com/bioinformatics-core-shared-training/Bulk_RNASeq_Course_March23)

## Extended materials

The [Extended Materials](Extended_index.md) contain extensions to some of the
sessions and additional materials, including instruction on downloading and
processing the raw data for this course, a link to an excellent R course, and
where to get further help after the course.

## Additional Resources

* [Bioconductor for relevant R packages](https://bioconductor.org/)
* [DESeq2 Vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)  
* [RNAseq Workflow](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)  
* [RStudio CheatSheets](https://rstudio.com/resources/cheatsheets/)

## Acknowledgements

This course is based on the course [RNAseq analysis in
R](http://combine-australia.github.io/2016-05-11-RNAseq/) prepared by [Combine
Australia](https://combine.org.au/) and delivered on May 11/12th 2016 in
Carlton. We are extremely grateful to the authors for making their materials
available; Maria Doyle, Belinda Phipson, Matt Ritchie, Anna Trigos, Harriet
Dashnow, Charity Law.

![](Bulk_RNAseq_Course_Base/images/combine_banner_small.png)

The materials have been rewritten/modified/corrected/updated by various
contributors over the past 5 years including:

Abigail Edwards
Ashley D Sawle
Chandra Chilamakuri
Dominique-Laurent Couturier
Guillermo Parada González
Hugo Tavares
Jon Price
Mark Dunning
Mark Fernandes
Oscar Rueda
Sankari Nagarajan
Stephane Ballereau
Tom Smith
Zeynep Kalender Atak

Apologies if we have missed anyone!

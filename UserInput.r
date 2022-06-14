
library(shiny)
library(dplyr)
library(data.table)



runApp(list(
    ui=fluidPage(
        textInput(inputId="covar",label="Covariate",placeholder="sex"),
        textInput(inputId="int_covar",label="Interaction Covariate"),
        radioButtons(inputId="int_type",label="Interaction Type",choiceNames=c('None','Addition','Multiplication'),choiceValues=c('','+','*')),
        radioButtons(inputId="term",label="Term",choiceNames=c('Linear','Product-Linear','Log-Linear'),choiceValues=c('lin','plin','loglin')),
        actionButton(inputId="submission",label="Add to list"),
        actionButton(inputId="ext",label="Exit Program"),
        textInput(inputId="file1",label="File Select",placeholder="temp.csv"),
        actionButton(inputId="load",label="Load Terms"),
        actionButton(inputId="save",label="Save Terms"),
        tableOutput("named_terms")
    ),
    server=function(input, output) {
        named_terms <- reactiveValues()
        named_terms$df <- data.table("item"=character(),"term"=character(),'covariate'=character(),'iteraction'=character(),'interaction covariate'=character())
        observeEvent(input$submission, {
            if (input$covar!=""){
                if (((input$int_type=='')&&(input$int_covar==''))||((input$int_type!='')&&(input$int_covar!=''))){
                    total_cov <- paste(input$covar,input$int_type,input$int_covar,sep="")
                    if (input$term=='lin'){
                        named_terms$df <- rbind(named_terms$df,list('item'=nrow(named_terms$df),'term'=paste("",total_cov,""),'covariate'=input$covar,'iteraction'=input$int_type,'interaction covariate'=input$int_covar))
                    } else if (input$term=='plin'){
                        named_terms$df <- rbind(named_terms$df,list('item'=nrow(named_terms$df),'term'=paste("1 + ",total_cov,""),'covariate'=input$covar,'iteraction'=input$int_type,'interaction covariate'=input$int_covar))
                    } else if (input$term=='loglin'){
                        named_terms$df <- rbind(named_terms$df,list('item'=nrow(named_terms$df),'term'=paste("e^(",total_cov,")"),'covariate'=input$covar,'iteraction'=input$int_type,'interaction covariate'=input$int_covar))
                    }
                }
                output$named_terms <- renderTable({named_terms$df})
            }
        })
        observeEvent(input$load, {
            if (file.exists(input$file1)){
                tb1 <- fread(file=input$file1,sep=",",blank.lines.skip=TRUE,header=TRUE,col.names=c("item","term","covariate","iteraction","interaction covariate"))
                tb1 <- tb1[item!='',]
                named_terms$df <- as.data.table(tb1)
            }
        })
        observeEvent(input$save, {
            if (input$file1!=""){
                if (grepl(".csv", input$file1)){
                    fwrite(named_terms$df,input$file1,sep=",",col.names=TRUE,row.names=FALSE)
                } else{
                    fwrite(named_terms$df,paste(input$file1,".csv",sep=""),sep=",",col.names=TRUE,row.names=FALSE,append=FALSE)
                }
            }
        })
        observeEvent(input$ext, {stopApp()})
    }))























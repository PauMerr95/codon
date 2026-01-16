#include <plog/Log.h>
#include <plog/Init.h>
#include <plog/Appenders/RollingFileAppender.h>
#include <plog/Formatters/CsvFormatter.h>

#ifndef PLOG_FILE_NAME
#error "Need to define PLOG_FILE_NAME"
#endif

static plog::RollingFileAppender<plog::CsvFormatter> fileAppender(PLOG_FILE_NAME);

struct LoggerInit {
    LoggerInit() {
        plog::init(plog::debug, &fileAppender);
    }
};


static LoggerInit loggerInit;

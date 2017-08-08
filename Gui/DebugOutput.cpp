/*********************************************
**
**	Created on: 	03/01/2014 2014
**	Author: 	matt - Matt Edmunds
**	File:		DebugOutput.cpp
**
**********************************************/

#include <iostream>
#include "DebugOutput.h"

DebugOutput* DebugOutput::sDebug = nullptr;

void myMessageOutput(QtMsgType type, const QMessageLogContext &/*context*/, const QString &msg)
{
	switch (type) {
    case QtInfoMsg:
        DebugOutput::Instance().printInfoMsg(msg);
        break;
	case QtDebugMsg:
		DebugOutput::Instance().printDebugMsg(msg);
		break;
	case QtWarningMsg:
		DebugOutput::Instance().printWarningMsg(msg);
		break;
	case QtCriticalMsg:
		DebugOutput::Instance().printCriticalMsg(msg);
		break;
	case QtFatalMsg:
		DebugOutput::Instance().printFatalMsg(msg);
		abort();
	}
}

DebugOutput& DebugOutput::Instance()
{
	if (sDebug == nullptr)
	{
		sDebug = new DebugOutput();
	}
	return *sDebug;
}

DebugOutput::DebugOutput()
{
	setupUi(this);

    qInstallMessageHandler(myMessageOutput);
}

DebugOutput::~DebugOutput()
{
    qInstallMessageHandler(0);

	if (sDebug != nullptr)
	{
		delete sDebug;
		sDebug = nullptr;
	}
}

void DebugOutput::printInfoMsg(const QString &msg)
{
    QString m("");
    m += msg;

    this->textBrowser->setTextColor( QColor( "black" ) );
    this->textBrowser->setFontItalic(false);
    this->textBrowser->append(m);

}

void DebugOutput::printDebugMsg(const QString &msg)
{
    QString m("Debug: ");
	m += msg;

	this->textBrowser->setTextColor( QColor( "black" ) );
	this->textBrowser->setFontItalic(false);
	this->textBrowser->append(m);

}

void DebugOutput::printWarningMsg(const QString &msg)
{
	QString m("Warning: ");
	m += msg;

	this->textBrowser->setTextColor( QColor( "black" ) );
	this->textBrowser->setFontItalic(true);
	this->textBrowser->append(m);

}

void DebugOutput::printCriticalMsg(const QString &msg)
{
	QString m("Critical: ");
	m += msg;

	this->textBrowser->setTextColor( QColor( "red" ) );
	this->textBrowser->setFontItalic(false);
	this->textBrowser->append(m);

}

void DebugOutput::printFatalMsg(const QString &msg)
{
	QString m("Fatal: ");
	m += msg;

	this->textBrowser->setTextColor( QColor( "red" ) );
	this->textBrowser->setFontItalic(true);
	this->textBrowser->append(m);
}

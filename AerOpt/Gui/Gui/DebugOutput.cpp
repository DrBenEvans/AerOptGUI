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

void myMessageOutput(QtMsgType type, const char *msg)
{
	switch (type) {
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

	qInstallMsgHandler(myMessageOutput);
}

DebugOutput::~DebugOutput()
{
	qInstallMsgHandler(0);

	if (sDebug != nullptr)
	{
		delete sDebug;
		sDebug = nullptr;
	}
}

void DebugOutput::printDebugMsg(const char *msg)
{
	QString m("Info: ");
	m += msg;

	this->textBrowser->setTextColor( QColor( "black" ) );
	this->textBrowser->setFontItalic(false);
	this->textBrowser->append(m);

}

void DebugOutput::printWarningMsg(const char *msg)
{
	QString m("Warning: ");
	m += msg;

	this->textBrowser->setTextColor( QColor( "black" ) );
	this->textBrowser->setFontItalic(true);
	this->textBrowser->append(m);

}

void DebugOutput::printCriticalMsg(const char *msg)
{
	QString m("Critical: ");
	m += msg;

	this->textBrowser->setTextColor( QColor( "red" ) );
	this->textBrowser->setFontItalic(false);
	this->textBrowser->append(m);

}

void DebugOutput::printFatalMsg(const char *msg)
{
	QString m("Fatal: ");
	m += msg;

	this->textBrowser->setTextColor( QColor( "red" ) );
	this->textBrowser->setFontItalic(true);
	this->textBrowser->append(m);
}

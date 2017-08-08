/*********************************************
**
**	Created on: 	01/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		Canvas.cpp
**
**********************************************/

#include <QDebug>
#include <QPainter>
#include <QMouseEvent>
#include <cmath>

#include "Canvas.h"
#include "ConstraintsDialog.h"
#include "Arrow.h"

Canvas::Canvas(QWidget *parent) :  QWidget(parent),
    mChordLength(1.0),
    mWidthMin(-0.75),
    mWidthMax(0.75),
    mHeightMin(-0.5),
    mHeightMax(0.5),
    mMesh(nullptr),
    mProfile(nullptr)
{
    colourmap.emplace_back( 0,         0,        0,        0.423499   );
	colourmap.emplace_back( 0.0668895, 0,        0.119341, 0.529244   );
	colourmap.emplace_back( 0.133779,  0,        0.238697, 0.634974   );
	colourmap.emplace_back( 0.200669,  0,        0.346853, 0.687877   );
	colourmap.emplace_back( 0.267559,  0,        0.450217, 0.718135   );
	colourmap.emplace_back( 0.334448,  0,        0.553552, 0.664836   );
	colourmap.emplace_back( 0.401338,  0,        0.651087, 0.51931    );
	colourmap.emplace_back( 0.468227,  0.115846, 0.724788, 0.35285    );
	colourmap.emplace_back( 0.535117,  0.326772, 0.781201, 0.140185   );
	colourmap.emplace_back( 0.602006,  0.522759, 0.79852,  0.0284581  );
	colourmap.emplace_back( 0.668897,  0.703166, 0.788678, 0.00885023 );
	colourmap.emplace_back( 0.735786,  0.845121, 0.751141, 0          );
	colourmap.emplace_back( 0.802675,  0.955734, 0.690822, 0          );
	colourmap.emplace_back( 0.869565,  0.995407, 0.56791,  0.0618448  );
	colourmap.emplace_back( 0.936455,  0.987716, 0.403403, 0.164858   );
	colourmap.emplace_back( 1,         0.980407, 0.247105, 0.262699   );

	setMouseTracking(true);
	installEventFilter(this);
	setContextMenuPolicy(Qt::CustomContextMenu);
    //Calculate initial canvas coordinates here.
	mHeightScale = 0;
	mWidthScale = 0;
	mCurrentMinHeight = 0;
	mCurrentMinWidth = 0;
	mOffset = mChordLength * 0.5;
	mHighlight = -1;
	calcCanvasScale();
}

void Canvas::setMesh(QSharedPointer<Mesh> mesh) {
    mMesh = mesh;
    update();
}

void Canvas::setProfile(QSharedPointer<Profile> profile) {
    mProfile = profile;
}

Canvas::~Canvas()
{
}

void Canvas::paintEvent(QPaintEvent *event)
{
	Q_UNUSED(event);
	QPainter painter(this);

	calcCanvasScale();
	drawBackground(painter);
    drawProfile(painter);
    drawMesh(painter);
	drawLogos(painter);
	drawScale(painter);
	drawAxis(painter);
}

bool Canvas::eventFilter(QObject* /*object*/, QEvent* event)
{
    //TODO - remove me
    return false;

	bool r = false;
	QPoint pos = {0,0};
	int loc = -1;

	if (event->type() == QEvent::MouseMove)
	{
		r = true;
		pos = static_cast<QMouseEvent*>(event)->pos();
        loc = pickNodeCheck(pos, mMesh);
		if (loc != mHighlight)
		{
			mHighlight = loc;
			repaint();
		}
	}

	if (event->type() == QEvent::MouseButtonDblClick)
	{
		r = true;
		pos = static_cast<QMouseEvent*>(event)->pos();
        loc = pickNodeCheck(pos, mMesh);
		if (loc != -1)
		{
			//store id here in data object
            mMesh->selectControlPoint(loc);
            BoundaryPoint& point = mMesh->getControlPoint(loc);
            point.setBoundCoords(0, 0, 0, 0);

			repaint();
		}
	}

	if (event->type() == QEvent::MouseButtonPress)
	{
		if( static_cast<QMouseEvent*>(event)->button() == Qt::RightButton )
		{
			pos = static_cast<QMouseEvent*>(event)->pos();
            loc = pickNodeCheck(pos, mMesh);
			if (loc != -1)
			{
				//find loc in getControlPoints() list
				//if exists then show menu
                auto& control = mMesh->getControlPoints();

                auto it = std::find(control.begin(), control.end(), loc);
                //TODO -> show menu item here...!

            }
		}
	}

	return r;
}

void Canvas::drawBackground(QPainter &painter)
{
	painter.setRenderHint(QPainter::Antialiasing, false);
	painter.setPen(palette().dark().color());
	painter.setBrush(Qt::white);
	painter.drawRect(QRect(0, 0, width() - 1, height() - 1));
}

void Canvas::drawLogos(QPainter& painter)
{
	QImage flag(":/images/AerOpt.png");
	painter.drawImage(QRect(1, 1, 105, 87), flag);

    //QImage coleg(":/images/ColegCenedlaetholCymraeg.png");
    //painter.drawImage(QRect(width()-150-1, height()-87-1, 150, 87), coleg);

    //QImage sulogo(":/images/SULogo.png");
    //painter.drawImage(QRect(width()-158-1, 1, 158, 100), sulogo);

	QImage map(":/images/Map.png");
	painter.drawImage(QRect(25, 100, 20, 300), map);

	painter.drawText( QPointF(47, 110), "Max");
	painter.drawText( QPointF(47, 400), "Min");
}

void Canvas::drawAxis(QPainter &painter)
{
	painter.setRenderHint(QPainter::Antialiasing, true);
	QPen pen(Qt::black, 1, Qt::SolidLine);
	painter.setBrush(QBrush(QColor(255, 255, 255, 0)));
	painter.setPen(pen);

    int arrow_length = 60;
    Arrow arrx( QPointF(30, height()-30), QPointF(arrow_length, height()-30) );
	arrx.paint(&painter);

    Arrow arry( QPointF(30, height()-30), QPointF(30, height()-arrow_length) );
	arry.paint(&painter);

    painter.drawText( QPointF(arrow_length+3, height()-30), "X");
    painter.drawText( QPointF(30, height()-(arrow_length+3)), "Y");
}

void Canvas::drawScale(QPainter &painter)
{
	painter.setRenderHint(QPainter::Antialiasing, true);
	QPen pen(Qt::black, 1, Qt::SolidLine);
	painter.setBrush(QBrush(QColor(255, 255, 255, 0)));
	painter.setPen(pen);

	painter.drawLine(w(0), h(-0.44), w(1), h(-0.44));
	painter.drawLine(w(0), h(-0.45), w(0), h(-0.42));
	painter.drawLine(w(0.5), h(-0.45), w(0.5), h(-0.43));
	painter.drawLine(w(1), h(-0.45), w(1), h(-0.42));

	painter.drawText( w(-0.007), h(-0.41), "0");
	painter.drawText( w(0.993), h(-0.41), "1");

}

void Canvas::drawProfile(QPainter &painter)
{
    if(!mProfile) {
        return;
    }

	painter.setRenderHint(QPainter::Antialiasing, false);
	QPen pen(Qt::black, 1, Qt::SolidLine);
	painter.setBrush(QBrush(QColor(255, 255, 255, 0)));
	painter.setPen(pen);
	QRect rect;

    for (auto& p : mProfile->getProfile())
	{
		rect.moveCenter( QPoint(w(p.first), h(p.second)) );
		rect.setWidth( 3 );
		rect.setHeight( 3 );
		painter.drawRect( rect );
	}
}

void Canvas::drawMesh(QPainter &painter)
{
    if(!mMesh) {
        return;
    }

    painter.setRenderHint(QPainter::Antialiasing, false);
	painter.setBrush(QBrush(QColor(255, 255, 255, 0)));
	painter.setPen(QPen(Qt::black, 1, Qt::SolidLine));
	QRect rect;
	QRect constraint;
	rect.setWidth( 7 );
	rect.setHeight( 7 );

    Boundaries& boundary = mMesh->getMeshBoundary();
    auto& bConnects = mMesh->getBConnects();
    auto& control = mMesh->getControlPoints();
    auto& meshpoints = mMesh->getMeshPoints();
    auto& meshconnects = mMesh->getMeshConnectivities();
    auto& resultsdata = mMesh->getMeshData();


	//This scope is for rendering the mesh when available
	painter.setPen( QPen(Qt::lightGray, 1, Qt::SolidLine) );
	painter.setRenderHint(QPainter::Antialiasing, true);
	if (meshconnects.size() > 0)
	{
		const uint r = 4;//pressure, 0 = rho, 1 = u, 2 = v, 3 = energy

		//Determin scalar range
		float min =  10000000;
		float max = -10000000;
		for (const auto& s : resultsdata)
		{
			const float& v = std::get<r>(s);
			if (v < min) min = v;
			if (v > max) max = v;
		}

		//Render triangulation
		for (auto& tuple : meshconnects)
		{
			const uint& i = std::get<0>(tuple);
			const uint& j = std::get<1>(tuple);
			const uint& k = std::get<2>(tuple);

			const std::pair<float,float>& p1 = meshpoints.at(i-1);
			const std::pair<float,float>& p2 = meshpoints.at(j-1);
			const std::pair<float,float>& p3 = meshpoints.at(k-1);

			if (resultsdata.size() == meshpoints.size())
			{
				const float& pressure1 = std::get<r>( resultsdata.at(i-1) );//Pressure
				const float& pressure2 = std::get<r>( resultsdata.at(j-1) );//Pressure
				const float& pressure3 = std::get<r>( resultsdata.at(k-1) );//Pressure
				const float pressure = (pressure1 + pressure2 + pressure3) * 0.333;

				const rgba&& colour = getColour(min, pressure, max);
				QColor qcolour(std::get<0>(colour),
							   std::get<1>(colour),
							   std::get<2>(colour),
							   std::get<3>(colour));

				QPainterPath path;
				path.moveTo( qreal(w( p1.first )), qreal(h( p1.second )) );
				path.lineTo( qreal(w( p2.first )), qreal(h( p2.second )) );
				path.lineTo( qreal(w( p3.first )), qreal(h( p3.second )) );
				path.lineTo( qreal(w( p1.first )), qreal(h( p1.second )) );
				painter.fillPath(path, QBrush(qcolour));
			}

			//Draw triangle boundary here
			painter.drawLine(
				w( p1.first ),
				h( p1.second ),
				w( p2.first ),
				h( p2.second )
				);

			painter.drawLine(
				w( p2.first ),
				h( p2.second ),
				w( p3.first ),
				h( p3.second )
				);

			painter.drawLine(
				w( p3.first ),
				h( p3.second ),
				w( p1.first ),
				h( p1.second )
				);
		}
	}

	//This scope is for rendering boundary points
    //from current boundary.
	painter.setRenderHint(QPainter::Antialiasing, false);
	painter.setPen( QPen(Qt::black, 1, Qt::SolidLine) );
    for (auto& p : boundary)
	{
        rect.moveCenter( QPoint( w(p.x()), h(p.y()) ) );
		rect.setWidth( 7 );
		rect.setHeight( 7 );
		painter.drawRect( rect );
	}

    //This scope is for rendering current and previous boundary
	painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setPen( QPen(Qt::blue, 2, Qt::SolidLine) );

    for (uint i = 0; i < bConnects.size(); ++i)
    {
        uint a = bConnects.at(i).first - 1;
        uint b = bConnects.at(i).second - 1;

        painter.drawLine(
                    w( boundary.at(a).x() ),
                    h( boundary.at(a).y() ),
                    w( boundary.at(b).x() ),
                    h( boundary.at(b).y() )
                    );
    }

    //This scope is for rendering boundary point connectivity
    //between all boundaries.
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setPen( QPen(Qt::green, 1, Qt::SolidLine) );
    if (mPrevMesh)
    {
        Boundaries& prevBoundaries = mPrevMesh->getMeshBoundary();
        const uint bSize = boundary.size();

        if (prevBoundaries.size() == bSize)
        {
            for (uint i = 0; i < bSize; i++)
            {
                prevBoundaries.at(i);
                boundary.at(i);

                painter.drawLine(
                            w( prevBoundaries.at(i).x() ),
                            h( prevBoundaries.at(i).y() ),
                            w( boundary.at(i).x() ),
                            h( boundary.at(i).y() )
                            );
            }
        }
    }

	//This scope is for control points only.
	painter.setRenderHint(QPainter::Antialiasing, false);
    painter.setPen( QPen(Qt::black, 1, Qt::SolidLine) );
    for (auto& i : control)
    {
        auto& p = boundary.at(i);

        //Set here red if selected, and green if fully constrained.
        if (p.isControlPoint())
        {
            painter.setBrush(QBrush(QColor(255, 0, 0, 255)));
        }
        else
        {
            painter.setBrush(QBrush(QColor(0, 255, 0, 255)));
        }

        rect.moveCenter( QPoint( w(p.x()), h(p.y()) ) );
        rect.setWidth( 7 );
        rect.setHeight( 7 );
        painter.drawRect( rect );
    }

    for (auto& i : control)
    {
        auto& p = boundary.at(i);
        if (!(p.isControlPoint()))
        {
            qreal x1,y1,x2,y2;
            p.getBoundCoords(&x1,&y1,&x2,&y2);

            constraint.setCoords( w(p.x()+x1), h(p.y()+y1), w(p.x()+x2), h(p.y()+y2) );

            painter.setBrush(QBrush(QColor(255, 0, 0, 69)));
            painter.drawRect( constraint );
        }
    }

    if (mHighlight > -1 && mHighlight < int(boundary.size()))
    {
        auto& p = boundary.at(mHighlight);
        //			painter.setBrush(QBrush(QColor(255, 215, 0, 180)));
        painter.setBrush(QBrush(QColor(255, 215, 0, 220)));
        rect.moveCenter( QPoint( w(p.x()), h(p.y()) ) );
        rect.setWidth( 7 );
        rect.setHeight( 7 );
        painter.drawRect( rect );
    }
}

int Canvas::h(float size) const
{
	float ret;

	size = size - currentMinHeight();
	ret = float( height() ) / heightScale();
	ret = ret * size;
	ret = std::round(ret);
	//Because height coordinates are inverted.
	ret = float( height() ) - ret;

	return ret;
}

int Canvas::w(float size) const
{
	float ret;

	size = size - currentMinWidth();
	size = size - offset();
	ret = float( width() ) / widthScale();
	ret = ret * size;
	ret = std::round(ret);

	return ret;
}

void Canvas::calcCanvasScale()
{
	float canvasRatio;
	float domainRatio;

	canvasRatio = float( width() ) / float( height() );
	domainRatio = ( widthMax() - widthMin() ) / ( heightMax() - heightMin() );

	if ( canvasRatio >= domainRatio )
	{
		//height is minimum
		mHeightScale = heightMax() - heightMin();
		mWidthScale = mHeightScale * canvasRatio;
		//now calc offsets to hold centre
		mCurrentMinHeight = heightMin();
		mCurrentMinWidth = widthMin() * ( mWidthScale / ( widthMax() - widthMin() ));
	}
	else
	{
		//width is minimum
		mWidthScale = widthMax() - widthMin();
		mHeightScale = mWidthScale * 1/canvasRatio;
		//now calc offsets to hold centre
		mCurrentMinWidth = widthMin();
		mCurrentMinHeight = heightMin() * ( mHeightScale / ( heightMax() - heightMin() ));
	}
}

int Canvas::pickNodeCheck(const QPoint& pos, QSharedPointer<Mesh> mesh)
{

    Boundaries& boundary = mesh->getMeshBoundary();
    if(boundary.size() > 0) {
        int index = -1;
        bool r = true;
        int i = 0;
        QRect rect;

        for (const auto& point : boundary)
        {
            r = true;

            rect.moveCenter( QPoint( w(point.x()), h(point.y()) ) );
            rect.setWidth( 13 );
            rect.setHeight( 13 );
            r &= rect.contains( pos );
            if (r) index = i;

            ++i;
        }
        return index;
    } else {
        return -1;
    }
}

const rgba Canvas::getColour(const float& min, const float& value, const float& max) const
{
	bool c = true;
	rgba&& colour = std::make_tuple(255, 255, 255, 255);
	float norm;

	c &= min <  max;

	if (c)
	{
		norm = (value - min) / (max - min);

		//compute rgba from norm
		transferFunction(norm, colour);
	}

	return colour;
}

void Canvas::transferFunction(const float& value, rgba& colour) const
{
	bool b = false;
	uint i;
	for (i = 0; i < colourmap.size()-1; i++)
	{
		b |= value < std::get<0>(colourmap.at(i+1));
		if (b) break;
	}

	if (i < colourmap.size()-1)
	{

		float diff = std::get<0>(colourmap.at(i+1)) - std::get<0>(colourmap.at(i+0));
		float v = value - std::get<0>(colourmap.at(i+0));
		float linear = v / diff;

		float diffr = std::get<1>(colourmap.at(i+1)) - std::get<1>(colourmap.at(i+0));
		float r = (diffr * linear) + std::get<1>(colourmap.at(i+0));

		float diffg = std::get<2>(colourmap.at(i+1)) - std::get<2>(colourmap.at(i+0));
		float g = (diffg * linear) + std::get<2>(colourmap.at(i+0));

		float diffb = std::get<3>(colourmap.at(i+1)) - std::get<3>(colourmap.at(i+0));
		float b = (diffb * linear) + std::get<3>(colourmap.at(i+0));

		std::get<0>(colour) = std::round(255.0 * r);
		std::get<1>(colour) = std::round(255.0 * g);
		std::get<2>(colour) = std::round(255.0 * b);
	}
}

void Canvas::setConstraints(const unsigned int index)
{
	//Open new dialogue here to get constraints.
    ConstraintsDialog diag(*mMesh, this);
	diag.setConstraint(index);
	diag.show();
	diag.exec();
}

void Canvas::resetConstraints(const unsigned int index)
{
    BoundaryPoint& point = mMesh->getControlPoint(index);
    point.setBoundCoords(0.0, 0.0, 0.0, 0.0);
}

//Getters and Setters
const float& Canvas::heightMax() const
{
	return mHeightMax;
}

const float& Canvas::heightMin() const
{
	return mHeightMin;
}

const float& Canvas::widthMax() const
{
	return mWidthMax;
}

const float& Canvas::widthMin() const
{
	return mWidthMin;
}

const float& Canvas::widthScale() const
{
	return mWidthScale;
}

const float& Canvas::heightScale() const
{
	return mHeightScale;
}

const float& Canvas::currentMinHeight() const
{
	return mCurrentMinHeight;
}

const float& Canvas::currentMinWidth() const
{
	return mCurrentMinWidth;
}

const float& Canvas::offset() const
{
	return mOffset;
}

#include "MeshGraphicsItem.h"
#include "BoundaryPointGraphicsItem.h"
#include <QtWidgets>

MeshGraphicsItem::MeshGraphicsItem(int scale, QGraphicsItem* parent) :
    mScale(scale)
{
    this->color = QColor(0,255,0);
    setZValue(-10);
    setPos(0,0);

    // set initial bounding box
    mXmax =  -std::numeric_limits<float>::infinity();
    mYmax =  -std::numeric_limits<float>::infinity();
    mXmin =  std::numeric_limits<float>::infinity();
    mYmin =  std::numeric_limits<float>::infinity();
}

void MeshGraphicsItem::setMesh(std::shared_ptr<Mesh> mesh) {
    mMesh = mesh;
    setBoundaryPoints(mMesh->getMeshBoundary());
    meshChanged();
}

void MeshGraphicsItem::meshChanged() {
    calcBoundingBox();
    prepareGeometryChange();
    QGraphicsItem::update();
}

void MeshGraphicsItem::showControlPoints(bool visible) {
    foreach(auto& item, this->childItems()) {
        item->setVisible(visible);
    }

}

void MeshGraphicsItem::calcBoundingBox() {
    if(!mMesh) return;

    float x, y;

    for (auto& p : mMesh->getMeshPoints())
    {
        x = p.first*mScale;
        y = p.second*mScale;

        if(x>mXmax)
            mXmax = x;
        if(x<mXmin)
            mXmin = x;
        if(y>mYmax)
            mYmax = y;
        if(y<mYmin)
            mYmin = y;

    }

    int margin=10;

    mBoundingBox = QRectF(mXmin-margin, mYmin-margin, (mXmax-mXmin)+2*margin, (mYmax-mYmin)+2*margin);
}

void MeshGraphicsItem::setBoundaryPoints(Boundaries& boundaryPoints) {
    foreach(auto& item, this->childItems()) {
        delete item;
    }

    // Draw control point objects, if mMesh is set
    for (auto& p : boundaryPoints)
    {
        BoundaryPointGraphicsItem* cp = new BoundaryPointGraphicsItem(mScale, this);
        cp->setPos(w(p.x()), h(p.y()));
    }
}

QRectF MeshGraphicsItem::boundingRect() const
{
    return mBoundingBox;
}

QPainterPath MeshGraphicsItem::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}

void MeshGraphicsItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    Q_UNUSED(widget);
    if(!mMesh) return;

    painter->setRenderHint(QPainter::Antialiasing, false);
    painter->setBrush(QBrush(QColor(255, 255, 255, 0)));
    painter->setPen(QPen(Qt::black, getBrushSize(), Qt::SolidLine));
    QRect rect;
    rect.setWidth( 7 );
    rect.setHeight( 7 );

    const auto& boundary = mMesh->getMeshBoundary();
    const auto& bConnects = mMesh->getBConnects();
    const auto& meshpoints = mMesh->getMeshPoints();
    const auto& meshconnects = mMesh->getMeshConnectivities();
    const auto& resultsdata = mMesh->getMeshData();


    //This scope is for rendering the mesh when available
    painter->setPen( QPen(Qt::lightGray, getBrushSize(), Qt::SolidLine) );
    painter->setRenderHint(QPainter::Antialiasing, true);
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

                QColor qcolour(255,0,0,0);

                QPainterPath path;
                path.moveTo( qreal(w( p1.first )), qreal(h( p1.second )) );
                path.lineTo( qreal(w( p2.first )), qreal(h( p2.second )) );
                path.lineTo( qreal(w( p3.first )), qreal(h( p3.second )) );
                path.lineTo( qreal(w( p1.first )), qreal(h( p1.second )) );
                painter->fillPath(path, QBrush(qcolour));
            }

            //Draw triangle boundary here
            painter->drawLine(
                        w( p1.first ),
                        h( p1.second ),
                        w( p2.first ),
                        h( p2.second )
                        );

            painter->drawLine(
                        w( p2.first ),
                        h( p2.second ),
                        w( p3.first ),
                        h( p3.second )
                        );

            painter->drawLine(
                        w( p3.first ),
                        h( p3.second ),
                        w( p1.first ),
                        h( p1.second )
                        );
        }
    }

    //This scope is for rendering initial and last boundaries
    painter->setRenderHint(QPainter::Antialiasing, true);
    uint inc = bConnects.size()-1;
    if (inc < 1) inc = 1;

    painter->setPen( QPen(Qt::red, getBrushSize()*2, Qt::SolidLine) );

    for (uint j = 0; j < bConnects.size(); ++j)
    {
        uint a = bConnects.at(j).first - 1;
        uint b = bConnects.at(j).second - 1;

        painter->drawLine(
                    w( boundary.at(a).x() ),
                    h( boundary.at(a).y() ),
                    w( boundary.at(b).x() ),
                    h( boundary.at(b).y() )
                    );
    }
}

uint MeshGraphicsItem::getBrushSize() {
    return 0;
}

qreal MeshGraphicsItem::w(qreal width) {
    return width*mScale;
}

qreal MeshGraphicsItem::h(qreal height) {
    return height*mScale;
}

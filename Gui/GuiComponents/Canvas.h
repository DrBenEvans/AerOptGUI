/*********************************************
**
**	Created on: 	01/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		Canvas.h
**
**********************************************/

#ifndef CANVAS_H
#define CANVAS_H

#include <QWidget>
#include <QEvent>
#include "Mesh.h"
#include "Profile.h"

//Abbreviate long type names
typedef std::tuple<int,int,int,int> rgba;

class OptimisationRun;
/**
 * @brief The Canvas class
 * This class is used for rendering
 * data to a window/canvas.
 */
class Canvas : public QWidget
{
	Q_OBJECT
public:
	/**
	 * @brief Canvas
     * @param data A reference to the OptimisationRun class.
     * This class depends on the OptimisationRun class.
	 */
    explicit Canvas(QWidget *dialog);
	~Canvas();
    void setProfile(QSharedPointer<Profile> profile);
    void setMesh(QSharedPointer<Mesh> mesh);

protected:
	/**
	 * @brief paintEvent
	 * @param event The paint event.
	 * Overide inherited paint event function from QWidget.
	 */
	void paintEvent(QPaintEvent *event);

private:
	/**
	 * @brief drawBackground
	 * @param painter A reference to the current painter context.
	 * This function draws the canvas background.
	 */
	void drawBackground(QPainter &painter);
	/**
	 * @brief drawLogos
	 * @param painter A reference to the current painter context.
	 * This function draws the canvas logos.
	 */
	void drawLogos(QPainter &painter);
	/**
	 * @brief drawAxis
	 * @param painter A reference to the current painter context.
	 * This function draws the canvas axis.
	 */
	void drawAxis(QPainter &painter);
	/**
	 * @brief drawScale
	 * @param painter A reference to the current painter context.
	 * This function draws the canvas scale.
	 */
	void drawScale(QPainter &painter);
	/**
	 * @brief drawProfile
	 * @param painter A reference to the current painter context.
	 */
    void drawProfile(QPainter &painter);
	/**
	 * @brief drawMesh
	 * @param painter A reference to the current painter context.
	 */
    void drawMesh(QPainter &painter);

	/**
	 * @brief widthMin
	 * @return Minimum width of the domain.
	 */
	const float& widthMin() const;
	/**
	 * @brief widthMax
	 * @return Maximum width of the domain.
	 */
	const float& widthMax() const;
	/**
	 * @brief heightMin
	 * @return Minmum height of the domain.
	 */
	const float& heightMin() const;
	/**
	 * @brief heightMax
	 * @return Maximum height of the domain.
	 */
	const float& heightMax() const;
	/**
	 * @brief heightScale
	 * @return Scale of domain to screan space height.
	 */
	const float& heightScale() const;
	/**
	 * @brief widthScale
	 * @return Scale of domain to screan space width.
	 */
	const float& widthScale() const;
	/**
	 * @brief currentMinWidth
	 * @return Current minmum width of the domain with
	 * reference to current aspect ratio.
	 */
	const float& currentMinWidth() const;
	/**
	 * @brief currentMinHeight
	 * @return Current minmum height of the domain with
	 * reference to current aspect ratio.
	 */
	const float& currentMinHeight() const;
	/**
	 * @brief offset
	 * @return Current offset from centre of profile.
	 */
	const float& offset() const;

	/**
	 * @brief calcCanvasScale
	 * Calculates the current canvas scales for width
	 * and height with reference to the domain to canvas
	 * scales and aspect ratio of the canvas.
	 */
	void calcCanvasScale();
	/**
	 * @brief pickNodeCheck
	 * @param pos The position of the mouse curser.
     * @param data Read only reference to the OptimisationRun class.
	 * @return The index of the node the mouse curser is over.
	 */
    int pickNodeCheck(const QPoint& pos, QSharedPointer<Mesh> mesh);
	/**
	 * @brief getColour
	 * @param min Range minimum.
	 * @param value Value to be interpolated within the range.
	 * @param max Range maximum.
	 * @return The rgba interpolated colour.
	 */
	const rgba getColour(const float& min, const float& value, const float& max) const;
	/**
	 * @brief transferFunction
	 * @param value The range normalised value to be interpolated.
	 * @param colour A refeence to the colour object to be set by the function.
	 */
	void transferFunction(const float& value, rgba& colour) const;

	/**
	 * @brief h
	 * @param size Global space y value.
	 * @return Screan space inverted y value.
	 * Converts global coordinates to screan space coordinates.
	 */
	int h(float size) const;
	/**
	 * @brief w
	 * @param size Global space x value.
	 * @returnScrean space x value.
	 * Converts global coordinates to screan space coordinates.
	 */
	int w(float size) const;

	const float mChordLength;
	const float mWidthMin;
	const float mWidthMax;
	const float mHeightMin;
	const float mHeightMax;
	float mHeightScale;
	float mWidthScale;
	float mCurrentMinWidth;
	float mCurrentMinHeight;
	float mOffset;
	int mHighlight;

	std::vector<std::tuple<float,float,float,float>> colourmap;

    QSharedPointer<Mesh> mMesh;
    QSharedPointer<Mesh> mPrevMesh;
    QSharedPointer<Profile> mProfile;

signals:

public slots:
	bool eventFilter(QObject* object, QEvent* event);
	void setConstraints(const unsigned int index);
	void resetConstraints(const unsigned int index);
};

#endif // CANVAS_H

#include "quad_tree_maps/MergeAlgorithm.h"

std::pair<MergeAlgorithm::transformation, MergeAlgorithm::transformation> MergeAlgorithm::updateGauss(const transformation &sample)
{
    transformation variance;
    transformation mean;

    if(samples.size() == M)
    {
        for(const auto &s : samples)
        {
            mean += s / M;
        }
        for(const auto &s : samples)
        {
            variance += (s - mean) * (s - mean);
        }
        variance = variance / M;
        rotate(samples.begin(), samples.begin() + 1, samples.end());
        samples.back() = sample;

        return std::pair<transformation, transformation>(mean, variance);
    }
    else
    {
        samples.push_back(sample);
        return std::pair<transformation, transformation>(transformation(), transformation());
    }
}

void MergeAlgorithm::randomAdaptiveWalk(int numSteps, const QuadMap &m1, QuadMap &m2)
{
    int k{0};
    int step{0};
    long double bestFitness = deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());
    long double currentFitness{0};
    long double bestbestFitness{INT32_MAX};
    double currentAI;
    double bestAI{0};
    double bestTransformAI{0};
    int bestIter;
    int itersSinceLastNewBest{0};
    int cellsInAgr;

    std::multimap<double, transformation, std::greater<>> goodTransformations;
    transformation bestTransformation;
    std::pair<transformation, transformation> noise;
    transformation t(0, 0, 0);
    transformation s;
    std::vector <transformation> emergencyTransformations;

    while(step < numSteps)
    {
        //  Tworzenie nowej próbki
        s.tX = t.tX + randGenerator(mean, std::max(minVariance, noise.second.tX));
        s.tY = t.tY + randGenerator(mean, std::max(minVariance, noise.second.tY));
        s.theta = t.theta + randGenerator(mean, std::max(minVariance, noise.second.theta));

        //  Transformacja mapy i obliczenie wartosci wskaznikow
        m2.transformMapInPlace(s.tX, s.tY, degToRad(s.theta));
        cellsInAgr = cellsInAgreement(*m1.getOriginalMap(), *m2.getMovingMap());
        currentFitness = deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());
        currentAI = overlap(*m1.getOriginalMap(), *m2.getMovingMap());

        if(currentAI == 0 || currentAI == 100 || cellsInAgr < 50)
        {
            t = emergencyTransformations.front();
            s.tX = t.tX + randGenerator(mean, std::max(minVariance, noise.second.tX));
            s.tY = t.tY + randGenerator(mean, std::max(minVariance, noise.second.tY));
            s.theta = t.theta + randGenerator(mean, std::max(minVariance, noise.second.theta));

            //  Transformacja mapy
            m2.transformMapInPlace(s.tX, s.tY, degToRad(s.theta));
            cellsInAgr = cellsInAgreement(*m1.getOriginalMap(), *m2.getMovingMap());
            currentFitness = deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());
            currentAI = overlap(*m1.getOriginalMap(), *m2.getMovingMap());
            samples.clear();
        }

        //  Akceptacja probki jezeli polepsza heurystyke, lub zostala wylosowana przez randomSelector()
        if(currentFitness < bestFitness || randomSelector(t, s) == s)
        {
            if(itersSinceLastNewBest == 1500) { std::cout << "ZAKONCZONO WCZESNIEJ\n"; break; }

            //  Jeżeli pokrycie jest odpowiednio wysokie, zapisz transformacje do pliku
            if(currentAI > 80 && currentAI < 99 && cellsInAgr > cellsTresh)
            {
                if(currentAI > bestAI) { bestAI = currentAI; }
                goodTransformations.insert(std::make_pair(currentAI * cellsInAgr, s));
            }

            // Jeżeli tylko nastąpiło polepszenie heurystyki
            if((currentFitness < bestbestFitness) && cellsInAgr > cellsTresh)
            {
                bestTransformation = s;
                bestbestFitness = currentFitness;
                bestTransformAI = currentAI;
                bestIter = k;
                itersSinceLastNewBest = 0;
            }
            ++k;
            ++itersSinceLastNewBest;
            t = s;
            bestFitness = currentFitness;

            //  Zapamietanie ostatnich probek
            if(emergencyTransformations.size() < 10)
            {
                emergencyTransformations.push_back(t);
            }
            else
            {
                rotate(emergencyTransformations.begin(), emergencyTransformations.begin() + 1, emergencyTransformations.end());
                emergencyTransformations.back() = t;
            }
            // Aktualizacja parametrow rozkladu Gaussa
            noise = updateGauss(s);
        }
        ++step;
    }
    //  Transformacja mapy zgodnie z najlepszym uzyskanym przeksztalceniem
    m2.transformMapInPlace(bestTransformation.tX, bestTransformation.tY, degToRad(bestTransformation.theta));
    std::cout << "\nNAJLEPSZE UZYSKANE PRZEKSZTALCENIE: " << bestTransformation << std::endl;
    std::cout << "NAJLEPSZA HEURYSTYKA: " << bestbestFitness << std::endl;
    std::cout << "WSKAZNIK NAJLEPSZEJ TRANSFORMACJI AI: " << bestTransformAI << std::endl;
    std::cout << "NAJLEPSZE REZULTATY UZYSKANO W ITERACJI: " << bestIter << std::endl;
    std::cout << "NAJLEPSZY WSKAZNIK AI: " << bestAI << std::endl;
    std::cout << "KONIEC\n";

    for(const auto &elem : goodTransformations)
    {
        std::cout << "AI * #CELLSAGREE = " << elem.first << " TRANSFORMATION: " << elem.second << std::endl;
    }
}

void MergeAlgorithm::newtonMinimization(const QuadMap &m1, QuadMap &m2, double tx, double ty, double theta)
{
    transformation t0 = transformation(tx, ty, theta);
    transformation firstOrder, secondOrder;
    transformation d, temp;
    double delta{0.2};

    for(int i = 0; i < 10; ++i)
    {
        firstOrder = transformation(0, 0, 0);
        secondOrder = transformation(0, 0, 0);
        d = transformation(0, 0, 0);

        m2.transformMapInPlace(t0.tX + delta, t0.tY, degToRad(t0.theta));
        firstOrder.tX += deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());
        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta));
        m2.transformMapInPlace(t0.tX - delta, t0.tY, degToRad(t0.theta));
        firstOrder.tX -= deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());

        m2.transformMapInPlace(t0.tX, t0.tY + delta, degToRad(t0.theta));
        firstOrder.tY += deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());
        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta));
        m2.transformMapInPlace(t0.tX, t0.tY - delta, degToRad(t0.theta));
        firstOrder.tY -= deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());

        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta + delta));
        firstOrder.theta += deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());
        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta));
        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta - delta));
        firstOrder.theta -= deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());

        firstOrder = firstOrder / (2 * delta);

        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta));
        m2.transformMapInPlace(t0.tX + delta, t0.tY, degToRad(t0.theta));
        secondOrder.tX += deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());
        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta));
        m2.transformMapInPlace(t0.tX - delta, t0.tY, degToRad(t0.theta));
        secondOrder.tX -= deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());
        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta));
        secondOrder.tX -= 2 * deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());

        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta));
        m2.transformMapInPlace(t0.tX, t0.tY + delta, degToRad(t0.theta));
        secondOrder.tY += deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());
        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta));
        m2.transformMapInPlace(t0.tX, t0.tY - delta, degToRad(t0.theta));
        secondOrder.tY -= deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());
        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta));
        secondOrder.tY -= 2 * deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());

        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta));
        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta + delta));
        secondOrder.theta += deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());
        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta));
        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta - delta));
        secondOrder.theta -= deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());
        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta));
        secondOrder.theta -= 2 * deltaHeuristic(*m1.getOriginalMap(), *m2.getMovingMap());

        secondOrder = secondOrder / (delta * delta);

        temp = firstOrder / secondOrder;
        d = -temp;

        t0 = t0 + d;
        m2.transformMapInPlace(t0.tX, t0.tY, degToRad(t0.theta));
    }
}

long double MergeAlgorithm::deltaHeuristic(const quadmap::QuadTree &m1, const quadmap::QuadTree &m2)
{
    return pictureDistance(m1, m2) - cLock * overlap(m1, m2);
}

MergeAlgorithm::transformation MergeAlgorithm::randomSelector(const transformation &t, const transformation &s)
{
    if(rand() % 3 == 0)
        return s;
    else
        return t;
}

double MergeAlgorithm::randGenerator(double mean, double variance)
{
    std::normal_distribution<double> distribution(mean, sqrt(variance));
    return distribution(gen);
}

std::vector<std::vector<int>> MergeAlgorithm::constructDMap(quadmap::QuadTree q, const cellState &state)
{
    minBoundingBox box = getBBX(q);
    std::vector<std::vector<int>> dmap(box.rows, std::vector<int>(box.columns, INT16_MAX));
    q.expand();
    int tempX{0};
    int tempY{0};
    int h{0};

    // Przypisuje wartosci 0 elementom w dMapie, ktore maja okreslony stan state
    for(auto it = q.begin_leafs(), end = q.end_leafs(); it!= end; ++it)
    {
        tempX = abs(box.topLeft.y() - it.getCoordinate().y()) / q.getResolution();
        tempY = (int)(box.columns - 1 - (abs(box.botRight.x() - it.getCoordinate().x())
                / q.getResolution()));

        if(binaryToState(q.isNodeOccupied(*it)) == state && tempX < dmap.size() && tempY < dmap.at(0).size())
        {
            dmap.at(tempX).at(tempY) = 0;
        }
    }
    for (int x = 1; x < box.rows; ++x) //  Start z lewego gornego naroznika
    {
        for (int y = 1; y < box.columns; ++y)
        {
            h = std::min(dmap.at(x - 1).at(y) + 1, dmap.at(x).at(y - 1) + 1);
            dmap.at(x).at(y) = std::min(dmap.at(x).at(y), h);
        }
    }
    for (int x = box.rows - 2; x >= 0; --x) //  Start z prawego dolnego naroznika
    {
        for (int y = box.columns - 2; y >= 0; --y)
        {
            h = std::min(dmap.at(x + 1).at(y) + 1, dmap.at(x).at(y + 1) + 1);
            dmap.at(x).at(y) = std::min(dmap.at(x).at(y), h);
        }
    }
    for (int x = box.rows - 2; x >= 0; --x) //  Start z lewego dolnego naroznika
    {
        for (int y = 1; y < box.columns; ++y)
        {
            h = std::min(dmap.at(x).at(y - 1) + 1, dmap.at(x + 1).at(y) + 1);
            dmap.at(x).at(y) = std::min(dmap.at(x).at(y), h);
        }
    }
    for (int x = 1; x < box.rows; ++x)  //  Start z prawego gornego naroznika
    {
        for (int y = box.columns - 2; y >= 0; --y)
        {
            h = std::min(dmap.at(x - 1).at(y) + 1, dmap.at(x).at(y + 1) + 1);
            dmap.at(x).at(y) = std::min(dmap.at(x).at(y), h);
        }
    }

    return dmap;
}

int MergeAlgorithm::dFunction(quadmap::QuadTree m1, const quadmap::QuadTree &m2, MergeAlgorithm::cellState state)
{
    int d{0};
    int tempX{0};
    int tempY{0};
    auto dmap = constructDMap(m2, state);
    minBoundingBox box = getBBX(m1);
    m1.expand();

    for(auto it = m1.begin_leafs(), end = m1.end_leafs(); it!= end; ++it)
    {
        if(binaryToState(m1.isNodeOccupied(*it)) == state)
        {
            tempX = abs(box.topLeft.y() - it.getCoordinate().y()) / m1.getResolution();
            tempY = (int)(box.columns - 1 - (abs(box.botRight.x() - it.getCoordinate().x())
                    / m1.getResolution()));
            if(tempX < dmap.size() && tempY < dmap.at(0).size())
            {
                d += dmap.at(tempX).at(tempY);
            }
        }
    }

    return d;
}

int MergeAlgorithm::pictureDistance(const quadmap::QuadTree &m1, const quadmap::QuadTree &m2)
{
    return dFunction(m2, m1, MergeAlgorithm::FREE) + dFunction(m2, m1, MergeAlgorithm::OCCUPIED)
        + dFunction(m1, m2, MergeAlgorithm::FREE) + dFunction(m1, m2, MergeAlgorithm::OCCUPIED);
}

MergeAlgorithm::minBoundingBox MergeAlgorithm::getBBX(const quadmap::QuadTree &q)
{
    quadmap::point2d tempTopLeft((float)INT32_MAX, (float)INT32_MIN);
    quadmap::point2d tempBotRight((float)INT32_MIN, (float)INT32_MAX);
    for(auto it = q.begin_leafs(), end = q.end_leafs(); it!= end; ++it)
    {
        if(it.getCoordinate().x() < tempTopLeft.x())
            tempTopLeft.x() = it.getCoordinate().x();

        else if(it.getCoordinate().x() > tempBotRight.x())
            tempBotRight.x() = it.getCoordinate().x();

        if(it.getCoordinate().y() > tempTopLeft.y())
            tempTopLeft.y() = it.getCoordinate().y();

        else if(it.getCoordinate().y() < tempBotRight.y())
            tempBotRight.y() = it.getCoordinate().y();
    }
    minBoundingBox box(tempTopLeft, tempBotRight);
    box.rows = abs(tempBotRight.y() - tempTopLeft.y()) / q.getResolution() + 1;
    box.columns = abs(tempBotRight.x() - tempTopLeft.x()) / q.getResolution() + 1;

    return box;
}

int MergeAlgorithm::cellsInAgreement(quadmap::QuadTree m1, const quadmap::QuadTree &m2)
{
    int agree{0};
    m1.expand();
    for(auto it = m1.begin_leafs(), end = m1.end_leafs(); it!= end; ++it)
    {
        quadmap::QuadTreeNode *n = m2.search(it.getCoordinate());
        if(n)
            if(binaryToState(m2.isNodeOccupied(n)) == binaryToState(m1.isNodeOccupied(*it)))
            {
                ++agree;
            }
    }

    return agree;
}

int MergeAlgorithm::cellsInDisagreement(quadmap::QuadTree m1, const quadmap::QuadTree &m2)
{
    int disagree{0};
    m1.expand();
    for(auto it = m1.begin_leafs(), end = m1.end_leafs(); it!= end; ++it)
    {
        quadmap::QuadTreeNode *n = m2.search(it.getCoordinate());
        if(n)
            if(binaryToState(m2.isNodeOccupied(n)) != binaryToState(m1.isNodeOccupied(*it)))
            {
                ++disagree;
            }
    }

    return disagree;
}

double MergeAlgorithm::overlap(const quadmap::QuadTree &m1, const quadmap::QuadTree &m2)
{
    if((double)cellsInAgreement(m1, m2) == 0)
        return 0.0;
    else
        return 100 * (double)cellsInAgreement(m1, m2)
        / (double)(cellsInAgreement(m1, m2) + cellsInDisagreement(m1, m2));
}

int MergeAlgorithm::numOfStates(const quadmap::QuadTree &q, bool isOccupied)
{
    int counter{0};
    for(auto it = q.begin_leafs(), end = q.end_leafs(); it != end; ++it)
    {
        if(q.isNodeOccupied(*it) == isOccupied)
        {
            ++counter;
        }
    }

    return counter;
}

MergeAlgorithm::cellState MergeAlgorithm::getCellState(const quadmap::QuadTree &q, const quadmap::point2d &p)
{
    quadmap::QuadTreeNode *n = q.search(p);
    if(!n)
        return UNKNOWN;
    else
        return binaryToState(q.isNodeOccupied(n));
}

double MergeAlgorithm::md(const quadmap::point2d &p1, const quadmap::point2d &p2)
{
    return abs(p1.x() - p2.x()) + abs(p1.y() - p2.y());
}

double MergeAlgorithm::md(const double &p1x, const double &p1y, const double &p2x, const double &p2y)
{
    return md(quadmap::point2d(p1x, p1y), quadmap::point2d(p2x, p2y));
}

MergeAlgorithm::cellState MergeAlgorithm::binaryToState(bool binary)
{
    return binary ? MergeAlgorithm::OCCUPIED : MergeAlgorithm::FREE;
}

bool MergeAlgorithm::stateToBinary(MergeAlgorithm::cellState state)
{
    return state == MergeAlgorithm::OCCUPIED;
}

double MergeAlgorithm::degToRad(const double &deg)
{
    return deg * M_PI / 180;
}





